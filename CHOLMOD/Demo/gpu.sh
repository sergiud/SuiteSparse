#!/bin/bash
CHOLMOD_USE_GPU=1 ../build/cholmod_l_demo < ~/nd6k.mtx
CHOLMOD_USE_GPU=0 ../build/cholmod_l_demo < ~/nd6k.mtx
