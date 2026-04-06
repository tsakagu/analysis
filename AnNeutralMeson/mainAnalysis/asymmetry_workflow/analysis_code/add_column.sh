#!/bin/bash

cat > fill.awk <<'AWK'
# Read extra.dat first: store columns 2 and 3 sequentially
FNR==NR { a[++m] = $2; b[m] = $3; next }

# Return position of the nth occurrence of character c in string s (1-based), or 0 if not found
function nthpos(s, c, n,   i, p) {
  p = 0
  for (i = 1; i <= n; i++) {
    p = index(substr(s, p + 1), c)
    if (p == 0) return 0
    p += (length(substr(s, 1, p + (p>0?0:0))) - length(substr(s, 1, p)))  # no-op; kept simple below
  }
  return 0
}

# Simpler nthpos implementation (no tricks)
function nthamp(s, n,   i, p, q) {
  q = 0
  for (i = 1; i <= n; i++) {
    p = index(substr(s, q + 1), "&")
    if (p == 0) return 0
    q += p
  }
  return q
}

{
  line = $0

  # Only touch lines that look like table rows: contain & and end with \\
  if (index(line, "&") && line ~ /\\\\[[:space:]]*$/) {
    row++

    if (row <= m) {
      # Find the 4th '&' (after col1..col4). This anchors "column 5 starts after here".
      p4 = nthamp(line, 4)
      if (p4 > 0) {
        prefix = substr(line, 1, p4)        # includes the 4th '&'
        tail   = substr(line, p4 + 1)       # starts at col5 content

        # Replace the first two empty cells in tail: " <spaces>&<spaces>&"
        # (these are the separators after empty col5 and empty col6)
        sub(/^[[:space:]]*&[[:space:]]*&/, " " a[row] " & " b[row] " ", tail)

        line = prefix tail
      }
    }
  }

  print line
}
AWK

awk -f fill.awk extra.dat table.tex > table_filled.tex
