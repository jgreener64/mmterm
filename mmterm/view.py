# View proteins and trajectories in the terminal

import os
import sys
import argparse
import termios
import fcntl
from io import StringIO
import logging

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
            # redirect stdin from pipe back to terminal
            sys.stdin = open('/dev/tty','r')
        except:
            logging.error("Piping structures not supported on this system (no /dev/tty)")
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
    info_str = "{} with {} models, {} chains ({}), {} residues, {} atoms".format(
                    os.path.basename(in_file), len(struc), len(chain_ids),
                    "".join(chain_ids), res_counter, atom_counter)

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

    # See https://stackoverflow.com/questions/13207678/whats-the-simplest-way-of-detecting-keyboard-input-in-python-from-the-terminal/13207724
    fd = sys.stdin.fileno()

    oldterm = termios.tcgetattr(fd)
    newattr = termios.tcgetattr(fd)
    newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
    termios.tcsetattr(fd, termios.TCSANOW, newattr)

    oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)

    canvas = Canvas()
    trans_x, trans_y = 0.0, 0.0
    rot_x, rot_y = 0.0, 0.0

    try:
        while True:
            os.system("clear")
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

            print(info_str)
            print("W/A/S/D rotates, T/F/G/H moves, I/O zooms, U spins, P cycles models, Q quits")
            canvas.clear()
            for x, y in points:
                canvas.set(x, y)
            print(canvas.frame())

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
                    k = sys.stdin.read(1)
                    if k:
                        if k.upper() == "O":
                            zoom /= zoom_speed
                        elif k.upper() == "I":
                            zoom *= zoom_speed
                        elif k.upper() == "F":
                            trans_x -= trans_speed
                        elif k.upper() == "H":
                            trans_x += trans_speed
                        elif k.upper() == "G":
                            trans_y -= trans_speed
                        elif k.upper() == "T":
                            trans_y += trans_speed
                        elif k.upper() == "S":
                            rot_x -= rot_speed
                        elif k.upper() == "W":
                            rot_x += rot_speed
                        elif k.upper() == "A":
                            rot_y -= rot_speed
                        elif k.upper() == "D":
                            rot_y += rot_speed
                        elif k.upper() == "U":
                            auto_spin = not auto_spin
                        elif k.upper() == "P" and len(struc) > 1:
                            cycle_models = not cycle_models
                        elif k.upper() == "Q":
                            return
                        break
                except IOError:
                    pass
    finally:
        termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)
        fcntl.fcntl(fd, fcntl.F_SETFL, oldflags)
