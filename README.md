# Molecule Viewer

<img src="MoleculeViewer.png" alt="Water molecule on a dark blue background" width="96">

Molecule Viewer is a desktop 3D viewer for molecular structures and trajectories in Extended XYZ and HDF5. It reads local files and loads coordinates on demand as you browse.

## Quick start

Use Python 3.12 in a virtual environment.

```sh
python -m pip install -r requirements.txt
python viewer.py water-demo.extxyz
```

[The demo file](water-demo.extxyz) contains two water structures with illustrative energy values. You can also start the viewer without arguments, then drag a file into the window or choose **Open file**.

## Using the viewer

- Select a molecular formula on the left and browse its conformers with the slider, **Previous** and **Next** buttons, or the conformer number field. Frames with the same elemental composition are grouped together, including isomers.
- Drag the molecule to rotate it and scroll to zoom. Toggle **Show hydrogen atoms** to include hydrogens in the view.
- Enable **Frame data** to see available energies, forces, charges, and other fields. Select a field to inspect its value or array components. Large arrays show a preview and report their total size.
- If HDF5 auto-detection is ambiguous, open **HDF5 fields** and select the coordinate and element datasets. You can use a glob such as `/molecules/*/coords` for repeated groups. Element and optional frame-offset paths can be relative to the coordinate group.

The viewer supports HDF5 coordinate arrays with shared or per-frame elements, packed arrays with frame offsets, and SPICE groups with atom-mapped SMILES. SPICE bonds come from the SMILES. Other bonds are inferred from distances. Coordinate units default to angstroms when the file does not specify them.

The first scan writes a local index. Input files are opened read-only, and the index is reused when an unchanged file is reopened.

## Building desktop packages

### macOS Apple Silicon

```sh
python -m pip install -r requirements.txt pyinstaller==6.22.3
python -m PyInstaller --noconfirm molecule_viewer.spec
```

The bundle is created at `dist/MoleculeViewer.app`.

### Windows x64

1. Extract the [Python 3.12.10 embeddable package](https://www.python.org/downloads/release/python-31210/) into a `runtime` directory next to the launcher. Copy `viewer.py`, `storage.py`, `hdf_layout.py`, `frame_data.py`, and `MoleculeViewer.png` into `runtime`.
2. From a regular Windows Python 3.12 installation, run `python -m pip install -r requirements.txt --target runtime/site-packages`. Copy `msvcp140*.dll` from `runtime/site-packages/PySide6` into `runtime`.
3. Set `runtime/python312._pth` to the following:

   ```text
   python312.zip
   .
   site-packages
   import site
   ```

4. With the Windows SDK and Zig 0.15.2 available, compile the icon resource and launcher:

   ```cmd
   rc.exe /fo viewer-icon.res windows_icon.rc
   zig cc --target=x86_64-windows-gnu -municode -Wl,--subsystem,windows -O2 -s windows_launcher.c viewer-icon.res -luser32 -o MoleculeViewer.exe
   ```

Keep `MoleculeViewer.exe` and the `runtime` folder together. The Windows launcher has been built, but execution on Windows has not been verified.

## License

The project source is available under the [MIT License](LICENSE). Dependencies retain their own licenses. The water-molecule icon was generated with ImageGen.
