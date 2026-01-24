# =============================================================================
# Render Professional HTML Report
# =============================================================================

cat("Generating professional HTML report...\n")

# Check for base64enc package
if (!requireNamespace("base64enc", quietly = TRUE)) {
  install.packages("base64enc", repos = "https://cloud.r-project.org")
}
library(base64enc)

# Read markdown content
md_lines <- readLines("HLA-I_PTM_Analysis_Report.md", warn = FALSE)

# Function to encode image to base64
encode_image <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw_data <- readBin(path, "raw", file.info(path)$size)
  base64 <- base64encode(raw_data)
  paste0("data:image/png;base64,", base64)
}

# Parse markdown and build HTML
html_parts <- list()
in_table <- FALSE
table_rows <- c()

for (i in seq_along(md_lines)) {
  line <- md_lines[i]

  # Skip empty lines
  if (trimws(line) == "") {
    if (in_table) {
      html_parts <- c(html_parts, paste0("<table>", paste(table_rows, collapse = ""), "</table>"))
      in_table <- FALSE
      table_rows <- c()
    }
    next
  }

  # Horizontal rule
  if (grepl("^---+$", line)) {
    html_parts <- c(html_parts, "<hr>")
    next
  }

  # Headers
  if (grepl("^# ", line)) {
    text <- sub("^# ", "", line)
    html_parts <- c(html_parts, paste0("<h1>", text, "</h1>"))
    next
  }
  if (grepl("^## ", line)) {
    text <- sub("^## ", "", line)
    html_parts <- c(html_parts, paste0("<h2>", text, "</h2>"))
    next
  }
  if (grepl("^### ", line)) {
    text <- sub("^### ", "", line)
    html_parts <- c(html_parts, paste0("<h3>", text, "</h3>"))
    next
  }

  # Images
  if (grepl("^!\\[", line)) {
    alt <- sub("^!\\[([^\\]]*)\\].*", "\\1", line)
    path <- sub(".*\\(([^)]+)\\).*", "\\1", line)
    base64_data <- encode_image(path)
    if (!is.null(base64_data)) {
      html_parts <- c(html_parts, paste0(
        '<figure><img src="', base64_data, '" alt="', alt, '">',
        '<figcaption>', alt, '</figcaption></figure>'
      ))
      cat("  Embedded:", path, "\n")
    }
    next
  }

  # Tables
  if (grepl("^\\|", line)) {
    if (!in_table) {
      in_table <- TRUE
      table_rows <- c()
    }
    # Skip separator row
    if (grepl("^\\|[-:|]+\\|$", line)) next

    cells <- strsplit(line, "\\|")[[1]]
    cells <- cells[cells != ""]
    cells <- trimws(cells)

    if (length(table_rows) == 0) {
      # Header row
      row <- paste0("<thead><tr>", paste0("<th>", cells, "</th>", collapse = ""), "</tr></thead><tbody>")
    } else {
      row <- paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
    }
    table_rows <- c(table_rows, row)
    next
  }

  # Close table if we hit non-table content
  if (in_table) {
    html_parts <- c(html_parts, paste0("<table>", paste(table_rows, collapse = ""), "</tbody></table>"))
    in_table <- FALSE
    table_rows <- c()
  }

  # Bold text
  line <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", line)

  # Regular paragraph
  html_parts <- c(html_parts, paste0("<p>", line, "</p>"))
}

# Close any remaining table
if (in_table) {
  html_parts <- c(html_parts, paste0("<table>", paste(table_rows, collapse = ""), "</tbody></table>"))
}

html_body <- paste(html_parts, collapse = "\n")

# Professional HTML template
html_doc <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HLA-I PTM Analysis Report</title>
    <style>
        @page { size: A4; margin: 2cm; }

        * { box-sizing: border-box; }

        body {
            font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
            font-size: 11pt;
            line-height: 1.5;
            color: #1a1a1a;
            max-width: 210mm;
            margin: 0 auto;
            padding: 40px;
            background: #fff;
        }

        h1 {
            font-size: 24pt;
            font-weight: 600;
            color: #1a365d;
            margin: 0 0 8px 0;
            padding-bottom: 12px;
            border-bottom: 3px solid #2c5282;
        }

        h1 + p {
            font-size: 12pt;
            color: #4a5568;
            margin-top: 0;
        }

        h2 {
            font-size: 16pt;
            font-weight: 600;
            color: #2d3748;
            margin: 32px 0 16px 0;
            padding-bottom: 8px;
            border-bottom: 2px solid #e2e8f0;
            page-break-after: avoid;
        }

        h3 {
            font-size: 12pt;
            font-weight: 600;
            color: #4a5568;
            margin: 24px 0 12px 0;
            page-break-after: avoid;
        }

        p {
            margin: 0 0 12px 0;
            text-align: justify;
        }

        hr {
            border: none;
            border-top: 1px solid #e2e8f0;
            margin: 24px 0;
        }

        table {
            width: 100%%;
            border-collapse: collapse;
            margin: 16px 0;
            font-size: 10pt;
        }

        th {
            background: #2c5282;
            color: white;
            font-weight: 600;
            text-align: left;
            padding: 10px 12px;
        }

        td {
            padding: 8px 12px;
            border-bottom: 1px solid #e2e8f0;
        }

        tr:nth-child(even) { background: #f7fafc; }

        figure {
            margin: 20px 0;
            text-align: center;
            page-break-inside: avoid;
        }

        img {
            max-width: 100%%;
            height: auto;
            border: 1px solid #e2e8f0;
            border-radius: 4px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }

        figcaption {
            font-size: 9pt;
            color: #718096;
            margin-top: 8px;
            font-style: italic;
        }

        strong { color: #2d3748; }

        @media print {
            body {
                padding: 0;
                font-size: 10pt;
            }
            h2 { page-break-before: auto; }
            figure { page-break-inside: avoid; }
            table { page-break-inside: avoid; }
        }
    </style>
</head>
<body>
%s
</body>
</html>', html_body)

# Write HTML file
writeLines(html_doc, "HLA-I_PTM_Analysis_Report.html")

cat("\nGenerated: HLA-I_PTM_Analysis_Report.html\n")
cat("File size:", round(file.info("HLA-I_PTM_Analysis_Report.html")$size / 1024 / 1024, 2), "MB\n")
cat("\nTo create PDF: Open HTML in browser and print to PDF\n")
