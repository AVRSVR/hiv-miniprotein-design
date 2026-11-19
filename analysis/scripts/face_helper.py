"""
A `face` command for PyMOL.

`orient` aligns a selection's principal axes. For a patch on a roughly spherical
protein surface that leaves the patch edge-on and anything bound to it hidden
behind the bulk of the molecule. `face` instead puts the camera on the ray from
the molecule's centre of mass through the patch's centre of mass, so the patch
is seen straight on.

    face <patch selection>, <selection to frame>, [buffer]
"""
import numpy as np
from pymol import cmd


def face(patch, around, buffer=3.0):
    com = np.array(cmd.centerofmass(around))
    pc = np.array(cmd.centerofmass(patch))
    z = pc - com
    n = float(np.linalg.norm(z))
    if n < 1e-6:
        cmd.orient(around)
        cmd.zoom(around, float(buffer))
        return
    z = z / n
    up = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(up, z))) > 0.9:
        up = np.array([0.0, 1.0, 0.0])
    x = np.cross(up, z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    v = list(cmd.get_view())
    v[0:9] = [float(t) for t in list(x) + list(y) + list(z)]
    cmd.set_view(v)
    cmd.zoom(around, float(buffer))


cmd.extend("face", face)
