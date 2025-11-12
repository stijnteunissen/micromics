dynamic_plot_width <- function(data, factor_cols) {
  max_len <- 0
  for (col in factor_cols) {
    if (col %in% colnames(data)) {
      # Zorg ervoor dat we de waarden als karakter behandelen
      col_lengths <- nchar(as.character(data[[col]]))
      max_len <- max(max_len, max(col_lengths, na.rm = TRUE))
    }
  }
  # Stel de breedte in op basis van de maximale tekenlengte
  if (max_len < 10) {
    return(8)
  } else if (max_len < 14) {
    return(12)
  } else if (max_len < 18) {
    return(14)
  } else {
    return(16)
  }
}
