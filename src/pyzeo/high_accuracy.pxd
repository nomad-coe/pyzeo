from libcpp.string cimport string
from pyzeo.netstorage cimport ATOM_NETWORK 

cdef extern from "../sphere_approx.h":
    cdef void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, 
            string AccSetting)
