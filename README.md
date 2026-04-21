# SK2CM

Convertisseur natif C++ pour exporter des releves Splatking / ARKit LiDAR vers un modele COLMAP. Le binaire cible Linux, macOS et Windows via CMake.

## Prerequis

- CMake 3.20+
- compilateur C++20
  - Linux: GCC 11+ ou Clang 14+
  - macOS: Apple Clang recent
  - Windows: Visual Studio 2022 ou clang-cl recent

## Build

Depuis la racine du depot :

```bash
cmake -S . -B build
cmake --build build --config Release
```

Executables produits :

- Linux / macOS: `build/sk2colmap`
- Windows multi-config: `build/Release/sk2colmap.exe`

## Usage

Exemple :

```bash
./build/sk2colmap --input ./LidarSeries --output ./output_final --validate-exported-model
```

Options utiles :

- `--align-to-colmap-model <Path>` pour aligner l'export sur le repere d'un modele COLMAP de reference
- `--snap-shared-cameras-to-reference` pour forcer exactement les poses des images communes sur ce modele de reference
- `--no-copy-images` pour ne pas recopier les JPEG dans le dossier de sortie
- `--validate-exported-model` pour relire les `.txt` et `.bin` exportes et comparer leurs comptes
- `--validation-images <N>` pour changer le nombre d'images utilisees pour la validation depth
- `--validation-stride <N>` pour changer le sous-echantillonnage des points LiDAR pendant la validation
- `--depth-bucket-size <N>` pour controler le sous-echantillonnage spatial lors de la projection LiDAR
- `--max-observations-per-track <N>` pour limiter le nombre d'observations gardees par point 3D
- `--min-capture-gap <N>` pour imposer un ecart minimum entre captures d'un meme track
- `--min-baseline-m <meters>` et `--min-angle-deg <degrees>` pour exiger une base geometrique minimale
- `--base-depth-tolerance-m <meters>` et `--relative-depth-tolerance <ratio>` pour ajuster le filtrage depth

## Sorties

Le convertisseur ecrit :

- `output/sparse/0/cameras.bin`
- `output/sparse/0/cameras.txt`
- `output/sparse/0/images.bin`
- `output/sparse/0/images.txt`
- `output/sparse/0/points3D.bin`
- `output/sparse/0/points3D.txt`
- `output/images/`
- `output/lidar_dense_world_points.ply`
- `output/lidar_sparse_world_points.ply`
- `output/model_validation.txt` quand `--validate-exported-model` est active
- `output/alignment_summary.txt` quand `--align-to-colmap-model` est utilise
- `output/camera_pose_summary.txt`
- `output/validation_summary.txt`
- `output/reprojection_summary.txt`

## Etat actuel

- Les poses et intrinsics sont convertis vers les formats COLMAP texte et binaire.
- Les JPEG portrait stockes en `1440x1920` sont automatiquement remappes depuis les intrinsics / depth maps paysage d'origine.
- Un `camera_id` est cree par image car les intrinsics varient legerement au cours de la capture.
- `lidar_dense_world_points.ply` correspond au nuage Splatking brut de session.
- `lidar_sparse_world_points.ply` et `points3D` sont construits a partir de vrais points LiDAR, en ne gardant que les observations image / depth coherentes.
- `--snap-shared-cameras-to-reference` privilegie l'alignement visuel des cameras communes sur un modele COLMAP de reference, au prix possible d'une degradation de la coherence reprojection des points LiDAR.
- `camera_pose_summary.txt` permet de verifier numeriquement que les centres camera exportes correspondent bien aux poses ARKit source.
- Un raffinement de track elimine les observations trop redondantes et exige une base geometrique minimale entre vues.

## Notes

- La nouvelle implementation portable vit dans `src/` avec un build CMake a la racine.
- Validation recente sur le jeu `LidarSeries_20260420_105850` : meme nombre de cameras, d'images et de points 3D que la version C#, avec une seule observation 2D supplementaire conservee cote C++ grace au calcul en double precision et au correctif de round-to-even sur le mapping depth.
