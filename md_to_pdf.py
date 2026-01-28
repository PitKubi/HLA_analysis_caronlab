#!/usr/bin/env python3
"""Convert markdown report to PDF via HTML + Chrome headless."""
import markdown2
import subprocess
import sys
import os

md_file = "HLA-I_PTM_Analysis_Report.md"
html_file = "HLA-I_PTM_Analysis_Report.html"
pdf_file = "HLA-I_PTM_Analysis_Report.pdf"

with open(md_file, "r") as f:
    md_content = f.read()

html_body = markdown2.markdown(md_content, extras=["tables", "fenced-code-blocks"])

html_full = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  @page {{
    size: A4;
    margin: 20mm 15mm;
  }}
  body {{
    font-family: 'Helvetica Neue', Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.5;
    color: #222;
    max-width: 100%;
  }}
  h1 {{
    font-size: 20pt;
    border-bottom: 2px solid #333;
    padding-bottom: 6px;
    margin-top: 0;
  }}
  h2 {{
    font-size: 15pt;
    color: #1a1a1a;
    border-bottom: 1px solid #ccc;
    padding-bottom: 4px;
    page-break-after: avoid;
  }}
  h3 {{
    font-size: 12pt;
    color: #333;
    page-break-after: avoid;
  }}
  table {{
    border-collapse: collapse;
    width: auto;
    margin: 12px 0;
    font-size: 10pt;
  }}
  th, td {{
    border: 1px solid #999;
    padding: 6px 12px;
    text-align: left;
  }}
  th {{
    background-color: #f0f0f0;
    font-weight: bold;
  }}
  tr:nth-child(even) {{
    background-color: #f9f9f9;
  }}
  img {{
    max-width: 100%;
    height: auto;
    display: block;
    margin: 10px auto;
    page-break-inside: avoid;
  }}
  /* Make donut plot use full width */
  img[alt="Figure 1B"] {{
    width: 100%;
  }}
  /* Make circos plots use full width */
  img[alt="Figure 4A"] {{
    width: 100%;
  }}
  img[alt="Figure 4A2"] {{
    width: 100%;
  }}
  hr {{
    border: none;
    border-top: 1px solid #ddd;
    margin: 20px 0;
  }}
  p {{
    margin: 8px 0;
  }}
  strong {{
    color: #111;
  }}
</style>
</head>
<body>
{html_body}
</body>
</html>
"""

with open(html_file, "w") as f:
    f.write(html_full)

result = subprocess.run([
    "google-chrome",
    "--headless",
    "--disable-gpu",
    "--no-sandbox",
    "--print-to-pdf=" + os.path.abspath(pdf_file),
    "--print-to-pdf-no-header",
    os.path.abspath(html_file)
], capture_output=True, text=True)

if result.returncode == 0:
    size = os.path.getsize(pdf_file) / (1024*1024)
    print(f"Created {pdf_file} ({size:.1f} MB)")
else:
    print("Error:", result.stderr)
    sys.exit(1)

# Clean up HTML
os.remove(html_file)
