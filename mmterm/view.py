# View proteins and trajectories in the terminal

import os
import sys
import argparse
import curses
from io import StringIO

import numpy as np
from drawille import Canvas, line

def view_protein(in_file, file_format=None, curr_model=1, chains=[], box_size=100.0):
    if box_size < 10.0 or box_size > 400.0:
        print("Box size must be between 10 and 400")
        return

    zoom_speed = 1.1
    trans_speed = 1.0
    rot_speed = 0.1
    spin_speed = 0.01
    action_count = 500
    auto_spin = False
    cycle_models = False

    # Infer file format from extension
    if file_format is None:
        file_format = os.path.basename(in_file).rsplit(".", 1)[-1]

    # Handle stdin
    if in_file == "-":
        contents = sys.stdin.read()
        struct_file = StringIO(contents)
        try:
            # Redirect stdin from pipe back to terminal
            sys.stdin = open("/dev/tty", "r")
        except:
            print("Piping structures not supported on this system (no /dev/tty)")
            return
    else:
        struct_file = in_file

    if file_format.lower() == "pdb":
        from Bio.PDB import PDBParser
        p = PDBParser()
        struc = p.get_structure("", struct_file)
    elif file_format.lower() in ("mmcif", "cif"):
        from Bio.PDB.MMCIFParser import MMCIFParser
        p = MMCIFParser()
        struc = p.get_structure("", struct_file)
    elif file_format.lower() == "mmtf":
        from Bio.PDB.mmtf import MMTFParser
        struc = MMTFParser.get_structure(struct_file)
    else:
        print("Unrecognised file format")
        return

    # Get backbone coordinates
    coords = []
    connections = []
    atom_counter, res_counter = 0, 0
    chain_ids = []
    for mi, model in enumerate(struc):
        model_coords = []
        for chain in model:
            chain_id = chain.get_id()
            if len(chains) > 0 and chain_id not in chains:
                continue
            if mi == 0:
                chain_ids.append(chain_id)
            for res in chain:
                if mi == 0:
                    res_counter += 1
                res_n = res.get_id()[1]
                for atom in res:
                    if mi == 0:
                        atom_counter += 1
                    if atom.get_name() in (
                                    "N", "CA", "C", # Protein
                                    "P", "O5'", "C5'", "C4'", "C3'", "O3'", # Nucleic acid
                                ):
                        if mi == 0 and len(model_coords) > 0:
                            # Determine if the atom is connected to the previous atom
                            connections.append(chain_id == last_chain_id and (res_n == (last_res_n + 1) or res_n == last_res_n))
                        model_coords.append(atom.get_coord())
                        last_chain_id, last_res_n = chain_id, res_n
        model_coords = np.array(model_coords)
        if mi == 0:
            if model_coords.shape[0] == 0:
                print("Nothing to show")
                return
            coords_mean = model_coords.mean(0)
        model_coords -= coords_mean # Center on origin of first model
        coords.append(model_coords)
    coords = np.array(coords)
    if curr_model > len(struc):
        print("Can't find that model")
        return

    # Build help strings
    info_str = "{} with {} models, {} chains ({}), {} residues, {} atoms".format(
                    os.path.basename(in_file), len(struc), len(chain_ids),
                    "".join(chain_ids), res_counter, atom_counter)
    help_str = "W/A/S/D rotates, T/F/G/H moves, I/O zooms, U spins, P cycles models, Q quits"

    # Make square bounding box of a set size and determine zoom
    x_min, x_max = float(coords[curr_model - 1, :, 0].min()), float(coords[curr_model - 1, :, 0].max())
    y_min, y_max = float(coords[curr_model - 1, :, 1].min()), float(coords[curr_model - 1, :, 1].max())
    x_diff, y_diff = x_max - x_min, y_max - y_min
    box_bound = float(np.max([x_diff, y_diff])) + 2.0
    zoom = box_size / box_bound
    x_min = zoom * (x_min - (box_bound - x_diff) / 2.0)
    x_max = zoom * (x_max + (box_bound - x_diff) / 2.0)
    y_min = zoom * (y_min - (box_bound - y_diff) / 2.0)
    y_max = zoom * (y_max + (box_bound - y_diff) / 2.0)

    # Set up curses screen
    # https://docs.python.org/3/howto/curses.html
    stdscr = curses.initscr()
    curses.noecho()
    curses.cbreak()
    curses.curs_set(False)
    stdscr.keypad(True)

    # Divide curses screen into windows
    window_info = stdscr.subwin(2, curses.COLS - 1,  # height, width
                                0, 0)                # begin_y, begin_x
    window_structure = stdscr.subwin(curses.LINES - 1 - 2, curses.COLS - 1,  # height, width
                                     2, 0)                                   # begin_y, begin_x

    # Print help strings (only need to do this once)
    window_info.addnstr(0, 0, info_str, window_info.getmaxyx()[1] - 1)
    window_info.addnstr(1, 0, help_str, window_info.getmaxyx()[1] - 1)
    window_info.refresh()

    canvas = Canvas()
    trans_x, trans_y = 0.0, 0.0
    rot_x, rot_y = 0.0, 0.0

    try:
        while True:
            points = []

            for x_start, y_start, x_end, y_end in (
                                (x_min, y_min, x_max, y_min),
                                (x_max, y_min, x_max, y_max),
                                (x_max, y_max, x_min, y_max),
                                (x_min, y_max, x_min, y_min),
                            ):
                for x, y in line(x_start, y_start, x_end, y_end):
                    points.append([x, y])

            rot_mat_x = np.array([
                                    [1.0,           0.0,            0.0],
                                    [0.0, np.cos(rot_x), -np.sin(rot_x)],
                                    [0.0, np.sin(rot_x),  np.cos(rot_x)],
                                ], dtype=np.float32)
            rot_mat_y = np.array([
                                    [ np.cos(rot_y), 0.0, np.sin(rot_y)],
                                    [           0.0, 1.0,           0.0],
                                    [-np.sin(rot_y), 0.0, np.cos(rot_y)],
                                ], dtype=np.float32)
            trans_coords = coords[curr_model - 1] + np.array([trans_x, trans_y, 0.0], dtype=np.float32)
            zoom_rot_coords = zoom * np.matmul(rot_mat_y, np.matmul(rot_mat_x, trans_coords.T)).T

            for i in range(coords.shape[1] - 1):
                if connections[i]:
                    x_start, x_end = float(zoom_rot_coords[i, 0]), float(zoom_rot_coords[i + 1, 0])
                    y_start, y_end = float(zoom_rot_coords[i, 1]), float(zoom_rot_coords[i + 1, 1])
                    if x_min < x_start < x_max and x_min < x_end < x_max and y_min < y_start < y_max and y_min < y_end < y_max:
                        for x, y in line(x_start, y_start, x_end, y_end):
                            points.append([x, y])

            # Update displayed structure
            canvas.clear()
            for x, y in points:
                canvas.set(x, y)
            window_structure.addstr(0, 0, canvas.frame())
            window_structure.refresh()

            counter = 0
            while True:
                if auto_spin or cycle_models:
                    counter += 1
                    if counter == action_count:
                        if auto_spin:
                            rot_y += spin_speed
                        if cycle_models:
                            curr_model += 1
                            if curr_model > len(struc):
                                curr_model = 1
                        break
                try:
                    c = stdscr.getch()
                    if c != curses.ERR:
                        if c in (ord("o"), ord("O")):
                            zoom /= zoom_speed
                        elif c in (ord("i"), ord("I")):
                            zoom *= zoom_speed
                        elif c in (ord("f"), ord("F")):
                            trans_x -= trans_speed
                        elif c in (ord("h"), ord("H")):
                            trans_x += trans_speed
                        elif c in (ord("g"), ord("G")):
                            trans_y -= trans_speed
                        elif c in (ord("t"), ord("T")):
                            trans_y += trans_speed
                        elif c in (ord("s"), ord("S")):
                            rot_x -= rot_speed
                        elif c in (ord("w"), ord("W")):
                            rot_x += rot_speed
                        elif c in (ord("a"), ord("A")):
                            rot_y -= rot_speed
                        elif c in (ord("d"), ord("D")):
                            rot_y += rot_speed
                        elif c in (ord("u"), ord("U")):
                            auto_spin = not auto_spin
                        elif c in (ord("p"), ord("P")) and len(struc) > 1:
                            cycle_models = not cycle_models
                        elif c in (ord("q"), ord("Q")):
                            return
                        break
                except IOError:
                    pass
    except KeyboardInterrupt:
        # If user presses Ctrl+C, pretend as if they pressed q.
        return
    finally:
        # Teardown curses interface
        curses.nocbreak()
        curses.echo()
        curses.curs_set(True)
        stdscr.keypad(False)
        curses.endwin()

        # Make sure last view stays on screen
        print(info_str)
        print(help_str)
        print(canvas.frame())
