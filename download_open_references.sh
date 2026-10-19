#!/bin/sh
# Downloads the references that are legally free to download (public domain or public reports).
# Run on your own computer: sh download_open_references.sh
mkdir -p references && cd references
curl -L -o Griffith_1921_rupture_and_flow.pdf "https://ia802702.us.archive.org/31/items/philtrans01457975/01457975.pdf"
curl -L -o Ramberg_Osgood_1943_NACA_TN902.pdf "https://ntrs.nasa.gov/api/citations/19930081614/downloads/19930081614.pdf"
curl -L -o Rice_Rosengren_1967_report_HRR.pdf "https://apps.dtic.mil/sti/pdfs/AD0659300.pdf"
curl -L -o Krueger_2002_VCCT_NASA_CR-2002-211628.pdf "https://apps.dtic.mil/sti/pdfs/ADA402137.pdf"
echo "Done. Paywalled references: see references.html for DOI links."
