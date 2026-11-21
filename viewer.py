from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from math import cos, sin
from pathlib import Path
import json
import os
import sqlite3
import sys
import tempfile
from threading import Event

import numpy as np
from PySide6.QtCore import (QAbstractListModel, QEvent, QModelIndex, QObject, QPointF,
                            QSignalBlocker, QStandardPaths, Qt, QTimer, Signal)
from PySide6.QtGui import (QAction, QColor, QFont, QFontDatabase, QIcon, QKeySequence, QPainter, QPen)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QFileDialog, QFormLayout, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QListView, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton,
                               QSizePolicy, QSlider, QSplitter, QTreeWidget, QTreeWidgetItem,
                               QVBoxLayout, QWidget)

from frame_data import MAX_FIELDS, read_frame_fields
from hdf_layout import HdfMapping, HdfMappingNeeded, field_choices
from storage import FrameCache, Index, Source, build_index, load_frame


COLORS = {element: QColor(value) for element, value in {
    "H": "#e1e8ef", "N": "#6299fc", "O": "#ff6e78", "S": "#f4ca55",
    "P": "#ffab5e", "F": "#92e298", "Cl": "#68cb80", "Br": "#cf8d72",
    "I": "#bb95ff",
}.items()}
CARBON = (QColor("#779bbb"), QColor("#3bcab5"))
OTHER = QColor("#baabd9")
BACKGROUND = QColor("#121d2c")
TRACK_MAX = 1_000_000


def number(value):
    return f"{value:,}"


class Canvas(QWidget):
    def __init__(self):
        super().__init__()
        self.frame = None
        self.hydrogen = False
        self.pointer = None
        self.setMinimumSize(300, 260)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.reset()

    def reset(self):
        self.pitch, self.yaw, self.zoom = 0.32, -0.38, 1.0
        self.update()

    def show_frame(self, frame):
        self.frame = frame
        self.update()

    def show_hydrogen(self, state):
        self.hydrogen = state
        self.update()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.pointer = event.position()
            event.accept()

    def mouseReleaseEvent(self, event):
        self.pointer = None

    def mouseMoveEvent(self, event):
        if self.pointer is not None:
            delta = event.position() - self.pointer
            self.pointer = event.position()
            self.pitch += delta.y() * 0.012
            self.yaw += delta.x() * 0.012
            self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y() or event.pixelDelta().y()
        self.zoom = max(0.2, min(8.0, self.zoom * 1.12 ** (delta / 120)))
        self.update()
        event.accept()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), BACKGROUND)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        frame = self.frame
        if frame is None:
            painter.setPen(QColor("#bdccdc"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "Drop an .hdf5, .h5, or .extxyz file here")
            return

        cp, sp, cy, sy = cos(self.pitch), sin(self.pitch), cos(self.yaw), sin(self.yaw)
        rotation = np.array([[cy, 0, sy], [sp * sy, cp, -sp * cy], [-cp * sy, sp, cp * cy]])
        xyz = frame.xyz - frame.xyz.mean(axis=0)
        projected = xyz @ rotation.T
        scale = self.zoom * 0.40 * min(self.width(), self.height())
        scale /= max(1, np.linalg.norm(xyz, axis=1).max())
        x = self.width() / 2 + projected[:, 0] * scale
        y = self.height() / 2 - projected[:, 1] * scale
        z = projected[:, 2]
        visible = [symbol != "H" or self.hydrogen for symbol in frame.elements]
        if not any(visible):
            visible = [True] * len(visible)
        colors = [
            CARBON[min(frame.fragments[i], 1) if frame.fragments else 0]
            if symbol == "C" else COLORS.get(symbol, OTHER)
            for i, symbol in enumerate(frame.elements)
        ]
        shapes = [(z[i], True, i, i) for i, shown in enumerate(visible) if shown]
        shapes.extend(((z[a] + z[b]) / 2, False, a, b)
                      for a, b in frame.bonds if visible[a] and visible[b])
        for _, atom, a, b in sorted(shapes):
            point = QPointF(x[a], y[a])
            if not atom:
                other = QPointF(x[b], y[b])
                middle = (point + other) / 2
                for index, end in ((a, point), (b, other)):
                    pen = QPen(colors[index], max(2, min(16, scale * 0.10)))
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    painter.setPen(pen)
                    painter.drawLine(end, middle)
                continue
            symbol = frame.elements[a]
            radius = max(2, min(32, scale * (0.12 if symbol == "H" else 0.21)))
            painter.setPen(QPen(BACKGROUND))
            painter.setBrush(colors[a])
            painter.drawEllipse(point, radius, radius)
            if symbol not in {"C", "H"} and radius >= 7:
                font = QFont(painter.font())
                font.setPixelSize(max(8, int(radius)))
                font.setBold(True)
                painter.setFont(font)
                point.setX(point.x() - painter.fontMetrics().horizontalAdvance(symbol) / 2)
                point.setY(point.y() + painter.fontMetrics().ascent() / 3)
                painter.drawText(point, symbol)


class FormulaList(QAbstractListModel):
    page_size = 128

    def __init__(self):
        super().__init__()
        self.index_db = None
        self.total = 0
        self.page = []
        self.page_start = -1

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else self.total

    def assign(self, index_db):
        self.beginResetModel()
        self.index_db = index_db
        self.total = 0
        self.page_start = -1
        self.page = []
        self.endResetModel()

    def refresh(self, size):
        old = self.total
        self.page_start = -1
        if size > old:
            self.beginInsertRows(QModelIndex(), old, size - 1)
            self.total = size
            self.endInsertRows()
        if old:
            self.dataChanged.emit(self.index(0), self.index(old - 1))

    def entry(self, row):
        if self.index_db is None or not 0 <= row < self.total:
            return None
        start = row // self.page_size * self.page_size
        if start != self.page_start:
            try:
                self.page = self.index_db.formulas(start + 1, self.page_size)
            except sqlite3.Error:
                return None
            self.page_start = start
        return self.page[row - start] if row - start < len(self.page) else None

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or role not in (Qt.ItemDataRole.DisplayRole,
                                               Qt.ItemDataRole.AccessibleTextRole):
            return None
        row = self.entry(index.row())
        if row:
            return f"{row[1]}   ·   {number(row[2])}"


class MappingDialog(QDialog):
    def __init__(self, source, parent, reason=""):
        super().__init__(parent)
        self.setWindowTitle("HDF5 fields")
        self.resize(720, 240)
        self.automatic = False
        layout = QVBoxLayout(self)
        hint = QLabel(reason or "Choose the coordinate and element datasets.")
        hint.setTextFormat(Qt.TextFormat.PlainText)
        hint.setWordWrap(True)
        layout.addWidget(hint)
        help_text = QLabel(
            "Select a path or enter one. For repeated groups, "
            "use * in the coordinate path, such as /molecules/*/xyz. "
            "Element and offset paths can be relative to the coordinate group.")
        help_text.setWordWrap(True)
        layout.addWidget(help_text)

        form = QFormLayout()
        layout.addLayout(form)
        try:
            coordinates, elements, offsets = field_choices(source.path)
        except (OSError, ValueError) as error:
            coordinates = elements = offsets = []
            label = QLabel(str(error))
            label.setTextFormat(Qt.TextFormat.PlainText)
            label.setWordWrap(True)
            layout.addWidget(label)
        mapping = source.mapping
        self.coordinates = self.field(form, "Coordinates", coordinates,
                                      mapping.coordinates if mapping else "", "/data/positions")
        self.elements = self.field(form, "Symbols or atomic numbers", elements,
                                   mapping.elements if mapping else "", "atomic_numbers")
        self.offsets = self.field(form, "Frame offsets", offsets,
                                  mapping.offsets if mapping else "", "For flat arrays")
        self.units = QComboBox()
        for label, value in (("From file, defaults to Å", ""), ("Angstroms (Å)", "angstrom"),
                             ("Nanometers (nm)", "nm"), ("Bohr (bohr)", "bohr")):
            self.units.addItem(label, value)
        self.units.setCurrentIndex(self.units.findData(mapping.units if mapping else ""))
        form.addRow("Length units", self.units)
        buttons = QDialogButtonBox()
        self.apply = buttons.addButton("Open", QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.addButton("Cancel", QDialogButtonBox.ButtonRole.RejectRole)
        if mapping:
            automatic = buttons.addButton("Auto-detect", QDialogButtonBox.ButtonRole.ActionRole)
            automatic.clicked.connect(self.choose_auto)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.coordinates.currentTextChanged.connect(self.validate)
        self.elements.currentTextChanged.connect(self.validate)
        layout.addWidget(buttons)
        self.validate()

    @staticmethod
    def field(form, label, choices, value, placeholder):
        box = QComboBox()
        box.setEditable(True)
        box.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        box.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        box.setMinimumContentsLength(24)
        for path, description in choices:
            box.addItem(path)
            box.setItemData(box.count() - 1, description, Qt.ItemDataRole.ToolTipRole)
        box.setCurrentIndex(-1)
        box.lineEdit().setPlaceholderText(placeholder)
        box.setEditText(value)
        form.addRow(label, box)
        return box

    def validate(self):
        self.apply.setEnabled(bool(self.coordinates.currentText().strip()
                                   and self.elements.currentText().strip()))

    def choose_auto(self):
        self.automatic = True
        self.accept()

    def mapping(self):
        path = self.coordinates.currentText().strip()
        return HdfMapping("/" + path.lstrip("/"), self.elements.currentText().strip(),
                          self.offsets.currentText().strip(), self.units.currentData())


class FrameDataPanel(QWidget):
    def __init__(self):
        super().__init__()
        self.preferred_field = ""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 2, 0, 0)
        header = QHBoxLayout()
        self.status = QLabel()
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        header.addWidget(self.status, 1)
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("Filter fields")
        self.filter.textChanged.connect(self.filter_rows)
        header.addWidget(self.filter)
        layout.addLayout(header)
        split = QSplitter()
        layout.addWidget(split, 1)
        self.table = QTreeWidget()
        self.table.setHeaderLabels(["Field", "Value", "Units"])
        self.table.setRootIsDecorated(False)
        self.table.setUniformRowHeights(True)
        self.table.header().setStretchLastSection(False)
        self.table.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.header().resizeSection(1, 95)
        self.table.header().resizeSection(2, 75)
        self.table.currentItemChanged.connect(self.show_value)
        split.addWidget(self.table)
        self.value = QPlainTextEdit()
        self.value.setAccessibleName("Frame field value")
        self.value.setReadOnly(True)
        self.value.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        split.addWidget(self.value)
        split.setSizes([500, 320])

    def set_fields(self, fields=(), message=""):
        self.table.clear()
        self.value.clear()
        if not message:
            if fields:
                message = f"Frame fields · {number(len(fields))}"
                if len(fields) >= MAX_FIELDS:
                    message += " (first 512)"
            else:
                message = "No additional data is available for this frame."
        self.status.setText(message)
        for field in fields:
            item = QTreeWidgetItem([field.name, field.summary, field.units])
            item.setToolTip(0, field.name)
            item.setData(0, Qt.ItemDataRole.UserRole, field.detail)
            self.table.addTopLevelItem(item)
        self.filter_rows(self.filter.text())

    def filter_rows(self, query):
        query = query.strip().casefold()
        first = preferred = None
        for i in range(self.table.topLevelItemCount()):
            item = self.table.topLevelItem(i)
            match = query in f"{item.text(0)} {item.text(2)}".casefold()
            item.setHidden(not match)
            if match and first is None:
                first = item
            if match and item.text(0) == self.preferred_field:
                preferred = item
        current = self.table.currentItem()
        if current is None or current.isHidden():
            self.table.setCurrentItem(preferred or first)

    def show_value(self, current, previous):
        if current is None:
            self.value.clear()
            return
        label = current.text(0)
        self.preferred_field = label
        if current.text(2):
            label += f" · {current.text(2)}"
        self.value.setPlainText(f"{label}\n\n{current.data(0, Qt.ItemDataRole.UserRole)}")


class WorkerSignals(QObject):
    indexed = Signal(int, object)
    scan_finished = Signal(int, object, bool, str, bool)
    frame_ready = Signal(int, object, object, str)
    fields_ready = Signal(int, object, str)


def scan_file(job, source, cancelled, signals):
    stage = None
    try:
        source.cache_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            name = source.cache_path.read_text(encoding="ascii").strip()
            if Path(name).name == name and name.startswith("scan-") and name.endswith(".sqlite3"):
                existing = source.cache_path.with_name(name)
                complete, _, _ = Index(existing).stats()
                if complete:
                    signals.scan_finished.emit(job, existing, False, "", False)
                    return
        except (OSError, sqlite3.Error, UnicodeError, TypeError):
            pass
        if cancelled.is_set():
            return
        with tempfile.NamedTemporaryFile(dir=source.cache_path.parent, prefix="scan-",
                                         suffix=".sqlite3", delete=False) as handle:
            stage = Path(handle.name)
        complete = build_index(source, stage, cancelled,
                               lambda: signals.indexed.emit(job, stage))
        if complete:
            signals.scan_finished.emit(job, stage, True, "", False)
            stage = None
    except Exception as error:
        if not cancelled.is_set():
            signals.scan_finished.emit(job, None, False, str(error), isinstance(error, HdfMappingNeeded))
    finally:
        if stage is not None:
            stage.unlink(missing_ok=True)


def read_frame_task(request, key, source, index, signals):
    try:
        frame = load_frame(source, index.location(*key))
        signals.frame_ready.emit(request, key, frame, "")
    except Exception as error:
        signals.frame_ready.emit(request, key, None, str(error))


def read_fields_task(request, key, source, index, cancelled, signals):
    try:
        if cancelled.is_set():
            return
        fields = read_frame_fields(source, index.location(*key), cancelled)
        if not cancelled.is_set():
            signals.fields_ready.emit(request, fields, "")
    except Exception as error:
        if not cancelled.is_set():
            signals.fields_ready.emit(request, (), str(error))


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.resize(1050, 700)
        self.setAcceptDrops(True)
        self.source = None
        self.index_db = None
        self.job = self.request = 0
        self.group = None
        self.count = self.position = 0
        self.cancelled = Event()
        self.temporary = None
        self.frame_cache = FrameCache()
        self.scan_pool = ThreadPoolExecutor(max_workers=2)
        self.frame_pool = ThreadPoolExecutor(max_workers=1)
        self.frame_future = None
        self.data_pool = ThreadPoolExecutor(max_workers=1)
        self.data_future = None
        self.data_cancelled = Event()
        self.shown_request = None
        self.signals = WorkerSignals()
        self.signals.indexed.connect(self.indexed)
        self.signals.scan_finished.connect(self.scan_finished)
        self.signals.frame_ready.connect(self.frame_ready)
        self.signals.fields_ready.connect(self.fields_ready)
        self.frame_timer = QTimer(self, singleShot=True, interval=28)
        self.frame_timer.timeout.connect(self.request_frame)
        self.cache_dir = Path(QStandardPaths.writableLocation(
            QStandardPaths.StandardLocation.CacheLocation)) / "indices"

        outer = QWidget()
        main = QVBoxLayout(outer)
        self.setCentralWidget(outer)
        header = QHBoxLayout()
        button = QPushButton("Open file…")
        button.clicked.connect(self.choose_file)
        header.addWidget(button)
        self.filename = QLabel("No file open")
        self.filename.setTextFormat(Qt.TextFormat.PlainText)
        self.filename.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        header.addWidget(self.filename, 1)
        self.progress = QLabel()
        header.addWidget(self.progress)
        self.mapping_button = QPushButton("HDF5 fields…")
        self.mapping_button.clicked.connect(lambda: self.choose_mapping())
        self.mapping_button.hide()
        header.addWidget(self.mapping_button)
        main.addLayout(header)
        action = QAction("Open file…", self)
        action.setShortcut(QKeySequence.StandardKey.Open)
        action.triggered.connect(self.choose_file)
        self.addAction(action)

        split = QSplitter()
        main.addWidget(split, 1)
        left = QWidget()
        left.setMinimumWidth(180)
        left_layout = QVBoxLayout(left)
        left_layout.addWidget(QLabel("Formulas"))
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find a formula")
        self.search.returnPressed.connect(self.find_formula)
        left_layout.addWidget(self.search)
        self.model = FormulaList()
        self.formulas = QListView()
        self.formulas.setUniformItemSizes(True)
        self.formulas.setLayoutMode(QListView.LayoutMode.Batched)
        self.formulas.setBatchSize(256)
        self.formulas.setModel(self.model)
        self.formulas.selectionModel().currentChanged.connect(self.select_group)
        left_layout.addWidget(self.formulas)
        split.addWidget(left)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        self.details = QLabel()
        self.details.setTextFormat(Qt.TextFormat.PlainText)
        self.details.setWordWrap(True)
        right_layout.addWidget(self.details)
        self.canvas = Canvas()
        self.view_split = QSplitter(Qt.Orientation.Vertical)
        self.view_split.addWidget(self.canvas)
        self.data_panel = FrameDataPanel()
        self.view_split.addWidget(self.data_panel)
        self.view_split.setChildrenCollapsible(False)
        self.view_split.setStretchFactor(0, 3)
        self.view_split.setStretchFactor(1, 2)
        right_layout.addWidget(self.view_split, 1)
        self.data_panel.hide()
        display = QHBoxLayout()
        hydrogen = QCheckBox("Show hydrogen atoms")
        hydrogen.toggled.connect(self.canvas.show_hydrogen)
        display.addWidget(hydrogen)
        self.data_toggle = QCheckBox("Frame data")
        self.data_toggle.toggled.connect(self.toggle_fields)
        display.addWidget(self.data_toggle)
        display.addStretch()
        reset = QPushButton("Reset view")
        reset.clicked.connect(self.canvas.reset)
        display.addWidget(reset)
        right_layout.addLayout(display)

        nav = QHBoxLayout()
        self.prev_button = QPushButton("Previous")
        self.prev_button.clicked.connect(lambda: self.set_position(self.position - 1))
        nav.addWidget(self.prev_button)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.valueChanged.connect(self.slide_to)
        self.slider.setAccessibleName("Conformer")
        nav.addWidget(self.slider, 1)
        self.next_button = QPushButton("Next")
        self.next_button.clicked.connect(lambda: self.set_position(self.position + 1))
        nav.addWidget(self.next_button)
        nav.addWidget(QLabel("Conformer"))
        self.frame_input = QLineEdit()
        self.frame_input.setAccessibleName("Conformer number")
        self.frame_input.setFixedWidth(85)
        self.frame_input.editingFinished.connect(self.jump)
        nav.addWidget(self.frame_input)
        self.frame_total = QLabel()
        nav.addWidget(self.frame_total)
        right_layout.addLayout(nav)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        split.setSizes([240, 800])
        self.update_nav()

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open structures", "",
                                             "Structures (*.hdf5 *.h5 *.extxyz);;All files (*)")
        if path:
            self.open_file(path)

    def dragEnterEvent(self, event):
        if any(url.isLocalFile() for url in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        paths = [url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()]
        if len(paths) == 1:
            event.acceptProposedAction()
            self.open_file(paths[0])
        elif paths:
            self.statusBar().showMessage("Drop one file at a time.", 5500)

    def choose_mapping(self, reason=""):
        if self.source is None or self.source.is_xyz:
            return
        dialog = MappingDialog(self.source, self, reason)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        mapping = None if dialog.automatic else dialog.mapping()
        self.open_file(self.source.path, mapping, remember=False)

    def open_file(self, path, mapping=None, remember=True):
        try:
            source = Source.open(path, self.cache_dir, mapping)
            if remember and mapping is None and not source.is_xyz:
                try:
                    config = json.loads(source.cache_path.with_suffix(".mapping").read_text(encoding="utf-8"))
                    source = Source.open(path, self.cache_dir, HdfMapping(**config))
                except (OSError, ValueError, TypeError):
                    pass
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Could not open file", str(error))
            return
        self.cancelled.set()
        self.cancelled = Event()
        self.job += 1
        self.invalidate_frame()
        if self.temporary:
            self.temporary.unlink(missing_ok=True)
            self.temporary = None
        self.source = source
        self.mapping_button.setVisible(not source.is_xyz)
        self.index_db = None
        self.group = None
        self.count = self.position = 0
        self.frame_cache = FrameCache()
        self.model.assign(None)
        self.search.clear()
        self.filename.setText(source.path.name)
        self.filename.setToolTip(str(source.path))
        self.statusBar().clearMessage()
        self.progress.setText("Indexing…")
        self.details.clear()
        self.canvas.show_frame(None)
        self.canvas.reset()
        self.update_nav()
        self.scan_pool.submit(scan_file, self.job, source, self.cancelled, self.signals)

    def indexed(self, job, path):
        if job == self.job and path.exists():
            self.refresh_index(path)

    def scan_finished(self, job, path, staging, error, needs_mapping):
        if job != self.job:
            if staging:
                path.unlink(missing_ok=True)
            return
        if error:
            self.invalidate_frame()
            self.group = self.index_db = None
            self.count = self.position = 0
            self.model.assign(None)
            self.update_nav()
            self.canvas.show_frame(None)
            self.details.clear()
            self.progress.setText("Read error")
            self.statusBar().showMessage(error)
            if needs_mapping:
                self.choose_mapping(error)
            else:
                QMessageBox.warning(self, "Could not read file", error)
            return
        if staging:
            pointer = path.with_suffix(".pointer")
            try:
                pointer.write_text(path.name, encoding="ascii")
                os.replace(pointer, self.source.cache_path)
            except OSError:
                self.temporary = path
                pointer.unlink(missing_ok=True)
        if not self.source.is_xyz:
            try:
                config = Source.index_path(self.source.path, self.cache_dir).with_suffix(".mapping")
                if self.source.mapping is None:
                    config.unlink(missing_ok=True)
                else:
                    config.write_text(json.dumps(asdict(self.source.mapping)), encoding="utf-8")
            except OSError:
                pass
        self.refresh_index(path)

    def refresh_index(self, path):
        index = Index(path)
        try:
            complete, total, groups = index.stats()
        except sqlite3.Error:
            return
        self.index_db = self.model.index_db = index
        self.model.refresh(groups)
        if self.group is None and groups:
            self.formulas.setCurrentIndex(self.model.index(0))
        elif self.group is not None:
            entry = self.model.entry(self.group - 1)
            if entry:
                self.count = entry[2]
                self.update_nav()
        summary = f"Frames: {number(total)} · Formulas: {number(groups)}"
        self.progress.setText(summary if complete else f"Indexing… {summary}")

    def find_formula(self):
        if self.index_db is None:
            return
        prefix = self.search.text().strip()
        if not prefix or not prefix.isalnum():
            return
        row = self.index_db.find(prefix)
        if row:
            index = self.model.index(row[0] - 1)
            self.formulas.setCurrentIndex(index)
            self.formulas.scrollTo(index)
        else:
            self.statusBar().showMessage("Formula not found in the indexed frames.", 4000)

    def select_group(self, index, previous):
        if not index.isValid():
            return
        entry = self.model.entry(index.row())
        if entry:
            self.group, _, self.count = entry
            self.position = 0
            self.canvas.reset()
            self.set_position(1)

    def slide_to(self, value):
        if self.count > 1:
            limit = min(TRACK_MAX, self.count - 1)
            self.set_position(1 + round(value * (self.count - 1) / limit))

    def jump(self):
        try:
            position = int(self.frame_input.text())
        except ValueError:
            position = self.position
        self.set_position(position)

    def update_nav(self):
        with QSignalBlocker(self.slider):
            limit = min(TRACK_MAX, max(0, self.count - 1))
            self.slider.setRange(0, limit)
            if limit:
                self.slider.setValue(round((self.position - 1) * limit / (self.count - 1)))
        self.slider.setEnabled(self.count > 1)
        self.prev_button.setEnabled(self.position > 1)
        self.next_button.setEnabled(self.position < self.count)
        self.frame_input.setEnabled(bool(self.count))
        position = str(self.position) if self.position else ""
        if not self.frame_input.hasFocus() or not self.frame_input.isModified():
            self.frame_input.setText(position)
        self.frame_total.setText(f"/ {number(self.count)}" if self.count else "")

    def invalidate_frame(self):
        self.request += 1
        self.shown_request = None
        self.cancel_fields()
        if self.data_toggle.isChecked():
            self.data_panel.set_fields(message="Waiting for frame…")
        self.frame_timer.stop()
        if self.frame_future:
            self.frame_future.cancel()

    def set_position(self, desired):
        if not self.count:
            return
        position = min(self.count, max(1, desired))
        if position != self.position:
            self.position = position
            self.invalidate_frame()
            self.details.setText("Loading frame…")
            self.frame_timer.start()
        self.frame_input.setModified(False)
        self.update_nav()

    def request_frame(self):
        key = self.group, self.position
        try:
            self.source.check_unchanged()
        except (OSError, ValueError) as error:
            self.frame_ready(self.request, key, None, str(error))
            return
        frame = self.frame_cache.get(key)
        if frame is not None:
            self.frame_ready(self.request, key, frame, "")
        else:
            self.frame_future = self.frame_pool.submit(
                read_frame_task, self.request, key, self.source, self.index_db, self.signals)

    def frame_ready(self, request, key, frame, error):
        if request != self.request:
            return
        self.canvas.show_frame(frame)
        if error:
            self.details.setText("Could not read frame. " + error)
            return
        self.shown_request = request
        self.frame_cache.add(key, frame)
        entry = self.model.entry(self.group - 1)
        details = f"{frame.description} · atoms: {number(len(frame.elements))}"
        if entry:
            details = f"{entry[1]} · {details}"
        if frame.inferred_bonds:
            details += " · distance-based bonds"
        self.details.setText(details)
        if self.data_toggle.isChecked():
            self.request_fields()

    def cancel_fields(self):
        self.data_cancelled.set()
        if self.data_future:
            self.data_future.cancel()

    def toggle_fields(self, visible):
        self.data_panel.setVisible(visible)
        if visible:
            height = self.view_split.height()
            self.view_split.setSizes([max(260, height - 210), 210])
            self.request_fields()
        else:
            self.cancel_fields()
            self.data_panel.set_fields()

    def request_fields(self):
        self.cancel_fields()
        if self.shown_request != self.request or self.index_db is None:
            self.data_panel.set_fields(message="Select a frame.")
            return
        self.data_panel.set_fields(message="Loading frame data…")
        self.data_cancelled = Event()
        self.data_future = self.data_pool.submit(
            read_fields_task, self.request, (self.group, self.position),
            self.source, self.index_db, self.data_cancelled, self.signals)

    def fields_ready(self, request, fields, error):
        if request == self.request and self.data_toggle.isChecked():
            message = "Could not read frame data. " + error if error else ""
            self.data_panel.set_fields(fields, message)

    def closeEvent(self, event):
        self.cancelled.set()
        self.job += 1
        self.invalidate_frame()
        self.scan_pool.shutdown(wait=True, cancel_futures=True)
        self.frame_pool.shutdown(wait=True, cancel_futures=True)
        self.data_pool.shutdown(wait=True, cancel_futures=True)
        QApplication.sendPostedEvents(self, QEvent.Type.MetaCall)
        if self.temporary:
            self.temporary.unlink(missing_ok=True)
        super().closeEvent(event)


class App(QApplication):
    window = None

    def event(self, event):
        if event.type() == QEvent.Type.FileOpen and event.file() and self.window:
            self.window.open_file(event.file())
            return True
        return super().event(event)


def main():
    app = App([sys.argv[0]])
    app.setOrganizationName("MoleculeViewer")
    app.setApplicationName("MoleculeViewer")
    app.setWindowIcon(QIcon(str(Path(__file__).with_name("MoleculeViewer.png"))))
    app.window = Window()
    app.window.show()
    if len(sys.argv) > 1:
        app.window.open_file(sys.argv[1])
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
