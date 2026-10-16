# Fracture Mechanics Companion

**Author:** Dr Régis Kubler, Arts et Métiers Institute of Technology. Developed with the assistance of Claude Opus 5.5 (Anthropic), October 2026. Code: MIT · Documentation: CC BY 4.0.

Interactive web pages for an MSc course on fracture mechanics. Every figure is computed live in the browser, from the elastic solutions up to a small finite element model of a cracked body.

Open `index.html` to start. It links to all the pages and shows an animated centre crack with its plastic zones and a J-integral contour.

## Pages

| File | Content |
|---|---|
| `index.html` | Landing page with the animation and links to everything |
| `crack-tip-stress-intensity.html` | Stress concentration K<sub>t</sub> at an elliptical hole and its sharp-crack limit, the 1/√r singularity, the Williams near-tip field, handbook K, K from FE data, the fracture check, small-scale yielding |
| `crack-tip-angular-functions.html` | Angular functions f<sub>ij</sub>(θ) for modes I, II and III, mixed-mode loading, T-stress |
| `crack-tip-j-integral.html` | The J-integral: definition, live contour integration, path independence, J = K²/E′, HRR field, T-stress and constraint, J from a bend test |
| `crack-tip-fe-fracture.html` | K and J in finite element software: quarter-point elements, a live spider-web mesh solution, extraction methods, collapsed-element variants, an Abaqus how-to and useful links |
| `notes/stress-intensity-notes.pdf` | LaTeX lecture notes with the equations and figures of the pages |

## Folder structure

```
.
├── index.html
├── crack-tip-stress-intensity.html
├── crack-tip-angular-functions.html
├── crack-tip-j-integral.html
├── crack-tip-fe-fracture.html
├── assets/            thumbnails used on the index page
├── notes/
│   ├── stress-intensity-notes.pdf
│   ├── stress-intensity-notes.tex
│   ├── make_figures.py
│   └── figures/       figures included by the .tex file
├── .nojekyll          serve files as they are on GitHub Pages
├── LICENSE            MIT (code)
├── LICENSE-docs.md    CC BY 4.0 (documentation)
└── README.md
```

## Publishing on GitHub Pages

1. Create a new repository and push the contents of this folder to the `main` branch, with `index.html` at the root.
2. On GitHub, open **Settings → Pages**.
3. Under **Build and deployment**, choose **Deploy from a branch**, then select `main` and `/ (root)`, and save.
4. After a minute the site is live at `https://<your-user>.github.io/<repository>/`.

```bash
git init
git add .
git commit -m "Fracture mechanics companion"
git branch -M main
git remote add origin https://github.com/<your-user>/<repository>.git
git push -u origin main
```

## Viewing locally

The pages are plain, self-contained HTML files with no build step and no dependencies. Double-click `index.html`, or serve the folder:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

The only external resource is Google Fonts; without a network connection the pages fall back to system fonts.

## Rebuilding the notes

```bash
cd notes
python3 make_figures.py          # needs numpy and matplotlib
pdflatex stress-intensity-notes.tex
pdflatex stress-intensity-notes.tex
```

## Licence

- **Code** (the JavaScript in the HTML pages and `notes/make_figures.py`): MIT, see [`LICENSE`](LICENSE).
- **Documentation** (page text, equations, figures and the lecture notes): CC BY 4.0, see [`LICENSE-docs.md`](LICENSE-docs.md).

## Author

**Dr Régis Kubler**, Arts et Métiers Institute of Technology.
Developed with the assistance of Claude Opus 5.5 (Anthropic), October 2026.

## Cite as

> Kubler, R. (2026). *Fracture Mechanics Companion*. Arts et Métiers Institute of Technology. Zenodo. DOI: https://doi.org/10.5281/zenodo.23421418

```bibtex
@software{kubler_2026_fracture,
  author    = {Kubler, Régis},
  title     = {Fracture Mechanics Companion},
  year      = {2026},
  month     = oct,
  publisher = {Zenodo},
  institution = {Arts et Métiers Institute of Technology},
  doi       = {https://doi.org/10.5281/zenodo.23421418}
}
```
