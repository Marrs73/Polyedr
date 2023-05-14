#!/usr/bin/env -S python3 -B

from time import time
from tk_drawer import TkDrawer
from polyedr import Polyedr


tk = TkDrawer()
try:
    for name in ["test_planes1", "test_planes2", "test_planes3",
                 "test_planes4", "test_box1", "test_box2"]:
        print("=============================================================")
        print(f"Начало работы с полиэдром '{name}'")
        start_time = time()
        Polyedr(f"data/{name}.geom").draw(tk)
        delta_time = time() - start_time
        print(f"Изображение полиэдра '{name}' заняло {delta_time} сек.")
        input("Hit 'Return' to continue -> ")
except (EOFError, KeyboardInterrupt):
    print("\nStop")
    tk.close()
