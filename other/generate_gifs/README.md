# GIFs generation with Materials visualization

## 1. Setup materials

1.1. Generate materials in Mat3ra JSON format with [Materials Designer](https://mat3ra-materials-designer.netlify.app/) or use a script.

1.2. Place materials JSON files in the `structures` directory.

1.3. Name them with a short name that will appear at bottom left of GIF.

## 2. Start wave.js

2.1. Run wave.js (from GitHub: https://github.com/Exabyte-io/wave.js) locally (default port 3002 -- used in notebooks).

Or keep URL of the web deployment in the notebook.

## 3. Record material GIFs

3.1. Run `record_gifs_with_wave.ipynb` to generate GIFs.

This notebook will generate GIFs for all JSON materials in the `structures` directory and save them on the top level of this repo (because we can't control the downloads directory for IDE).

3.2. Wait until the GIFs are downloaded. Move selected ones to the `input` folder or allow the next notebook to move them automatically.

## 4. Add overlays and generate final GIFs

4.1. Store any media files (e.g. images, videos) you want to overlay on the GIFs in the `assets` directory.

4.2. Run `gif_processing.ipynb` to add overlays and generate final GIFs.

This notebook will move the GIFs from the top level to the `output` directory, removing any duplications (judging by the file name), and add overlays with the material names.
