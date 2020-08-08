# View proteins and trajectories in the terminal

import os
import sys
import argparse
import curses
from io import StringIO

import numpy as np
from drawille import Canvas, line

zoom_speed = 1.1   # Scale factor / keypress
trans_speed = 1.0  # Motion / keypress
rot_speed = 0.1    # Rad / keypress
spin_speed = 0.01  # Rad / frame

protein_bb_atoms = ["N", "CA", "C"]
nucleic_acid_atoms = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
atoms_of_interest = protein_bb_atoms + nucleic_acid_atoms

def get_coords_schrodinger(struc, chains):
    coords, info, connections  = [], {}, []
    atom_counter, res_counter = 0, 0
    chain_ids = []
    # Make sure that all models have the same number of atoms
    if len(set([st.atom_total for st in struc])) > 1:
        print("Multiple models with varying number of atoms "
              "are not supported.")
        return None, None

    for mi, model in enumerate(struc):
        model_coords = []
        for chain in model.chain:
            chain_id = chain.name
            if chains and chain_ids not in chains:
                continue
            if mi == 0:
                chain_ids.append(chain_id)
            prev_res = None
            for res in structure.get_residues_by_connectivity(chain):
                res_counter += 1 if mi == 0 else 0
                res_n = res.resnum
                for atom in res.atom:
                    atom_counter += 1 if mi == 0 else 0
                    if atom.pdbname.strip() in atoms_of_interest:
                        if mi == 0 and len(model_coords) > 0:
                            # Determine if the atom is connected to the previous atom
                            connections.append(chain_id == last_chain_id and (res_n == (last_res_n + 1) or res_n == last_res_n))
                        model_coords.append(atom.xyz)
                        last_chain_id, last_res_n = chain_id, res_n
        model_coords = np.array(model_coords)
        if mi == 0:
            if model_coords.shape[0] == 0:
                return None, None
            coords_mean = model_coords.mean(0)
        model_coords -= coords_mean # Center on origin of first model
        coords.append(model_coords)

    info["chain_ids"] = chain_ids
    info["model_coords"] = model_coords
    info["atom_counter"] = atom_counter
    info["res_counter"] = res_counter
    info["num_struc"] = len(struc)
    info["connections"] =  connections
    return coords, info

def get_coords_biopython(struc, chains):
    coords, info, connections = [], {}, []
    atom_counter, res_counter = 0, 0
    chain_ids = []
    for mi, model in enumerate(struc):
        model_coords = []
        for chain in model:
            chain_id = chain.get_id()
            if chains and chain_id not in chains:
                continue
            if mi == 0:
                chain_ids.append(chain_id)
            for res in chain:
                res_counter += 1 if mi == 0 else 0
                res_n = res.get_id()[1]
                for atom in res:
                    atom_counter += 1 if mi == 0 else 0
                    if atom.get_name() in atoms_of_interest:
                        if mi == 0 and len(model_coords) > 0:
                            # Determine if the atom is connected to the previous atom
                            connections.append(chain_id == last_chain_id and (res_n == (last_res_n + 1) or res_n == last_res_n))
                        model_coords.append(atom.get_coord())
                        last_chain_id, last_res_n = chain_id, res_n
        model_coords = np.array(model_coords)
        if mi == 0:
            if model_coords.shape[0] == 0:
                return None, None
            coords_mean = model_coords.mean(0)
        model_coords -= coords_mean # Center on origin of first model
        coords.append(model_coords)

    info["chain_ids"] = chain_ids
    info["model_coords"] = model_coords
    info["atom_counter"] = atom_counter
    info["res_counter"] = res_counter
    info["num_struc"] = len(struc)
    info["connections"] =  connections
    return coords, info

def read_inputs(in_file, file_format, curr_model, chains):
    # Infer file format from extension
    file_format = file_format or os.path.basename(in_file).rsplit(".", 1)[-1]

    # Handle stdin
    if in_file == "-":
        contents = sys.stdin.read()
        struct_file = StringIO(contents)
        try:
            # Redirect stdin from pipe back to terminal
            sys.stdin = open("/dev/tty", "r")
        except:
            print("Piping structures not supported on this system (no /dev/tty)")
            return None, None
    else:
        struct_file = in_file

    # Use Biopython parser by default
    get_coords = get_coords_biopython

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
    elif file_format.lower() in ("mae", "maegz"):
        from schrodinger import structure
        struc = list(structure.StructureReader(struct_file))
        get_coords = get_coords_schrodinger
    else:
        print("Unrecognised file format")
        return None, None

    coords, info = get_coords(struc, chains)

    if coords is None or curr_model > len(coords):
        print("Nothing to show")
        return None, None

    return np.array(coords), info

def view(in_file, file_format=None, curr_model=1, chains=[], box_size=100.0):
    if box_size < 10.0 or box_size > 400.0:
        print("Box size must be between 10 and 400")
        return

    auto_spin = False
    cycle_models = False

    coords, info = read_inputs(in_file, file_format, curr_model, chains)
    if coords is None:
        return

    # Build help strings
    info_str = (f"{os.path.basename(in_file)} with {info['num_struc']} models, "
                f"{len(info['chain_ids'])} chains ({''.join(info['chain_ids'])}) "
                f"{info['res_counter']} residues, {info['atom_counter']} atoms.")
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
    stdscr.keypad(True)  # Respond to keypresses w/o Enter
    stdscr.nodelay(True)  # Don't block while waiting for keypress

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
        points = []
        do_update = True
        while True:
            curses.napms(50)  # Delay a short while

            # Re-draw structure if needed
            if do_update:
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
                    if not info['connections'][i]:
                        continue
                    x_start, x_end = float(zoom_rot_coords[i, 0]), float(zoom_rot_coords[i + 1, 0])
                    y_start, y_end = float(zoom_rot_coords[i, 1]), float(zoom_rot_coords[i + 1, 1])
                    # Check if the bond fits in the box
                    if x_min < x_start < x_max and x_min < x_end < x_max and y_min < y_start < y_max and y_min < y_end < y_max:
                        for x, y in line(x_start, y_start, x_end, y_end):
                            points.append([x, y])

                # Update displayed structure
                canvas.clear()
                for x, y in points:
                    canvas.set(x, y)
                window_structure.addstr(0, 0, canvas.frame())
                window_structure.refresh()
                do_update = False

            # Prepare rotation/model selection for next time
            if auto_spin:
                rot_y += spin_speed
                do_update = True
            if cycle_models:
                curr_model += 1
                if curr_model > len(coords):
                    curr_model = 1
                do_update = True

            # Handle keypresses
            try:
                c = stdscr.getch()
                if c != curses.ERR:
                    do_update = True
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
                    elif c in (ord("p"), ord("P")) and len(coords) > 1:
                        cycle_models = not cycle_models
                    elif c in (ord("q"), ord("Q")):
                        return
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
