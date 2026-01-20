read_table_auto <- function(path,
                            sheet = 1,
                            sep = NULL,
                            row_names = NULL,     # NULL=不设；1=第一列；或列名/列号
                            check.names = FALSE,
                            stringsAsFactors = FALSE,
                            ...) {
  stopifnot(is.character(path), length(path) == 1, file.exists(path))

  ext <- tolower(tools::file_ext(path))

  # 读取成 data.frame
  df <- switch(
    ext,
    "csv"  = utils::read.csv(path, header = TRUE, check.names = check.names,
                             stringsAsFactors = stringsAsFactors, ...),

    "tsv"  = utils::read.delim(path, header = TRUE, sep = "\t", check.names = check.names,
                               stringsAsFactors = stringsAsFactors, ...),

    "txt"  = {
      # 如果用户没指定 sep：优先猜测 tab，再退回空白分隔
      if (is.null(sep)) {
        # 尝试 tab
        df1 <- try(utils::read.delim(path, header = TRUE, sep = "\t",
                                     check.names = check.names,
                                     stringsAsFactors = stringsAsFactors, ...),
                   silent = TRUE)
        if (!inherits(df1, "try-error") && ncol(df1) > 1) {
          df1
        } else {
          utils::read.table(path, header = TRUE, sep = "", check.names = check.names,
                            stringsAsFactors = stringsAsFactors, ...)
        }
      } else {
        utils::read.table(path, header = TRUE, sep = sep, check.names = check.names,
                          stringsAsFactors = stringsAsFactors, ...)
      }
    },

    "xlsx" =,
    "xls"  = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("To read Excel files, please install the 'readxl' package.")
      }
      readxl::read_excel(path, sheet = sheet, ...)
      # read_excel 返回 tibble；下面会 as.data.frame
    },

    stop("Unsupported file type: .", ext,
         " (supported: csv, tsv, txt, xlsx/xls)")
  )

  df <- as.data.frame(df, check.names = check.names, stringsAsFactors = stringsAsFactors)

  # 设置行名
  if (!is.null(row_names)) {
    if (is.numeric(row_names) && length(row_names) == 1) {
      rn <- df[[row_names]]
      df[[row_names]] <- NULL
      rownames(df) <- make.unique(as.character(rn))
    } else if (is.character(row_names) && length(row_names) == 1) {
      if (!row_names %in% colnames(df)) {
        stop("row_names='", row_names, "' not found in columns.")
      }
      rn <- df[[row_names]]
      df[[row_names]] <- NULL
      rownames(df) <- make.unique(as.character(rn))
    } else {
      stop("row_names must be NULL, a single column index, or a single column name.")
    }
  }

  df
}


read_counts_and_meta <- function(counts_path,
                                 meta_path,
                                 counts_row_names = 1,
                                 meta_sample_col = 1,
                                 meta_group_col  = 2,
                                 counts_sheet = 1,
                                 meta_sheet = 1,
                                 ...) {

  counts_df <- read_table_auto(
    counts_path,
    sheet = counts_sheet,
    row_names = counts_row_names,
    ...
  )

  meta_df <- read_table_auto(
    meta_path,
    sheet = meta_sheet,
    row_names = NULL,
    ...
  )

  # 规范 meta 列
  if (ncol(meta_df) < max(meta_sample_col, meta_group_col)) {
    stop("meta.data must have at least two columns (Sample, Group).")
  }
  meta <- data.frame(
    Sample = as.character(meta_df[[meta_sample_col]]),
    Group  = as.character(meta_df[[meta_group_col]]),
    stringsAsFactors = FALSE
  )

  # 去掉空值
  meta <- meta[!is.na(meta$Sample) & meta$Sample != "", , drop = FALSE]

  # counts 转 matrix（数值化）
  counts_mat <- as.matrix(counts_df)
  storage.mode(counts_mat) <- "double"

  if (is.null(colnames(counts_mat))) {
    stop("counts table must have column names (sample IDs).")
  }

  # 基础一致性检查（不强行重排）
  miss_in_meta <- setdiff(colnames(counts_mat), meta$Sample)
  miss_in_mat  <- setdiff(meta$Sample, colnames(counts_mat))

  list(
    counts = counts_mat,
    meta = meta,
    check = list(
      missing_in_meta = miss_in_meta,
      missing_in_counts = miss_in_mat,
      setequal = setequal(colnames(counts_mat), meta$Sample),
      identical_order = identical(colnames(counts_mat), meta$Sample)
    )
  )
}
