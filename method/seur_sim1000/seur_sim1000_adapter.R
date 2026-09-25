ar <- commandArgs(trailingOnly = TRUE)
get_arg <- function(x) {
  i <- match(x, ar)
  if (is.na(i) || i == length(ar)) stop("missing command-line value for ", x)
  ar[[i + 1L]]
}

root <- normalizePath(get_arg("--rep-root"), winslash = "/", mustWork = TRUE)
ad <- file.path(root, "adapter")
if (dir.exists(ad)) stop("refusing to overwrite adapter output: ", ad)
dir.create(file.path(ad, "multi"), recursive = TRUE)
dir.create(file.path(ad, "susiex"))
dat <- readRDS(file.path(root, "input", "input.rds"))
an <- c("EUR", "EAS", "AFR")
if (!identical(names(dat$R_list), an) || !identical(names(dat$z_list), an) ||
    !identical(as.numeric(dat$n), c(1000000, 150000, 150000)) ||
    dat$p != 500L || dat$L != 5L) stop("adapter input contract failed")
vn <- sprintf("v%04d", seq_len(dat$p))
variants <- data.frame(CHR = 22L, SNP = vn, BP = 50000000L + seq_len(dat$p),
                       A1 = "G", A2 = "A")
write.table(variants, file.path(ad, "variants.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(field = c("p", "L", "ancestry_order", "n"),
                       value = c(dat$p, dat$L, paste(an, collapse = ","),
                                 paste(dat$n, collapse = ","))),
            file.path(ad, "multi", "meta.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
rt <- list()
i <- 0L
for (a in an) {
  R <- dat$R_list[[a]]
  z <- dat$z_list[[a]]
  rp <- file.path(ad, "multi", paste0("R_", a, ".f64le"))
  zp <- file.path(ad, "multi", paste0("z_", a, ".f64le"))
  writeBin(as.numeric(R), rp, size = 8L, endian = "little")
  writeBin(as.numeric(z), zp, size = 8L, endian = "little")
  Rr <- matrix(readBin(rp, "double", n = dat$p^2, size = 8L,
                       endian = "little"), dat$p, dat$p)
  zr <- readBin(zp, "double", n = dat$p, size = 8L, endian = "little")
  if (!identical(Rr, unname(R)) || !identical(zr, unname(as.numeric(z)))) {
    stop("MultiSuSiE binary roundtrip failed: ", a)
  }
  i <- i + 1L
  rt[[i]] <- data.frame(ancestry = a, target = "MultiSuSiE",
                         object = "R", max_abs_error = max(abs(Rr - R)))
  i <- i + 1L
  rt[[i]] <- data.frame(ancestry = a, target = "MultiSuSiE",
                         object = "z", max_abs_error = max(abs(zr - z)))

  se <- 1 / sqrt(dat$n[[a]])
  beta <- z / sqrt(dat$n[[a]])
  pv <- 2 * stats::pnorm(-abs(z))
  sp <- file.path(ad, "susiex", paste0("sst_", a, ".txt"))
  hdr <- paste(c("CHR", "SNP", "BP", "A1", "A2", "BETA", "SE", "STAT", "P"),
               collapse = "\t")
  rows <- paste(variants$CHR, variants$SNP, variants$BP, variants$A1,
                variants$A2, sprintf("%.17g", beta), sprintf("%.17g", se),
                sprintf("%.17g", z), sprintf("%.17g", pv), sep = "\t")
  writeLines(c(hdr, rows), sp, useBytes = TRUE)
  ss <- read.delim(sp, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(ss) != dat$p || !identical(ss$SNP, vn) ||
      max(abs(ss$STAT - z)) > 1e-14) stop("SuSiEx summary roundtrip failed: ", a)
  lp <- file.path(ad, "susiex", paste0("LD_", a, ".ld.bin"))
  writeBin(as.numeric(t(R)), lp, size = 4L, endian = "little")
  R4 <- matrix(readBin(lp, "numeric", n = dat$p^2, size = 4L,
                       endian = "little"), dat$p, dat$p, byrow = TRUE)
  if (max(abs(R4 - R)) > 1e-6 || max(abs(R4 - t(R4))) > 1e-6) {
    stop("SuSiEx LD roundtrip failed: ", a)
  }
  write.table(data.frame(22L, vn, 0L, variants$BP, variants$A1, variants$A2),
              file.path(ad, "susiex", paste0("LD_", a, "_ref.bim")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(data.frame(chr = 22L, snp_id = vn, a1 = "G", a2 = "A", maf = .25),
              file.path(ad, "susiex", paste0("LD_", a, "_frq.frq")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  i <- i + 1L
  rt[[i]] <- data.frame(ancestry = a, target = "SuSiEx",
                         object = "summary", max_abs_error = max(abs(ss$STAT - z)))
  i <- i + 1L
  rt[[i]] <- data.frame(ancestry = a, target = "SuSiEx",
                         object = "R", max_abs_error = max(abs(R4 - R)))
}
write.table(do.call(rbind, rt), file.path(root, "status", "adapter_roundtrip.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
fs <- list.files(ad, recursive = TRUE, full.names = TRUE)
fs <- fs[!file.info(fs)$isdir]
rel <- gsub("\\\\", "/", substring(fs, nchar(ad) + 2L))
o <- order(rel)
fs <- fs[o]
rel <- rel[o]
ip <- file.path(root, "input", "input.rds")
mf <- data.frame(
  case = dat$case, case_id = dat$case_id, rep = dat$rep, seed = dat$seed,
  input_md5 = unname(tools::md5sum(ip)), p = dat$p, L = dat$L,
  ancestry_order = paste(an, collapse = ","), n = paste(dat$n, collapse = ","),
  path = rel,
  format = ifelse(grepl("f64le$", rel), "float64_column_major_or_vector",
           ifelse(grepl("ld.bin$", rel), "float32_row_major", "tab_text")),
  md5 = as.character(tools::md5sum(fs)),
  size = format(file.info(fs)$size, scientific = FALSE, trim = TRUE))
write.table(mf, file.path(ad, "adapter_manifest.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(status = "success", input_md5 = mf$input_md5[[1L]],
                       files = nrow(mf),
                       manifest_md5 = unname(tools::md5sum(file.path(ad, "adapter_manifest.tsv")))),
            file.path(root, "status", "adapter_status.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
