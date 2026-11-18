from pathlib import Path

assets = Path(SPECPATH)
a = Analysis(
    ['viewer.py'],
    pathex=[SPECPATH],
    datas=[(str(assets / 'MoleculeViewer.png'), '.')],
    excludes=['IPython', 'matplotlib', 'pytest', 'tkinter', 'PyQt5', 'PyQt6', 'PySide2'],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name='MoleculeViewer',
    console=False,
    upx=False,
)
dist = COLLECT(exe, a.binaries, a.datas, name='MoleculeViewer', upx=False)
app = BUNDLE(
    dist,
    name='MoleculeViewer.app',
    version='1.3.1',
    icon=str(assets / 'MoleculeViewer.icns'),
    bundle_identifier='org.moleculeviewer.desktop',
    info_plist={
        'CFBundleDisplayName': 'Molecule Viewer',
        'NSHighResolutionCapable': True,
        'CFBundleDocumentTypes': [
            {
                'CFBundleTypeName': 'Molecular structures',
                'CFBundleTypeRole': 'Viewer',
                'LSHandlerRank': 'Alternate',
                'CFBundleTypeExtensions': ['hdf5', 'h5', 'extxyz'],
            },
        ],
    },
)
