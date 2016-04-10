# This is the Makefile for the mex-file version of Chaco, to be dynamically 
# linked with Matlab using the interface in the "meshpart" toolbox.  
# See file "README" for more information.
#
# You may need to change the following variables:
# 
# MATLAB      Matlab's top-level directory, often "/usr/matlab" or similar.
# CHACO       Chaco's source directory, normally the "code" subdirectory
#             of the Chaco distribution.
# MLCHACO     Directory with Chaco files modified for Matlab, normally
#             the directory containing this Makefile.
# DEST_DIR    Directory for the executable mlchaco.mex, normally the 
#             "meshpart" directory, which is normally the parent of the 
#             directory containing this Makefile.
# CC          Choice of C compiler.  gcc seems to work best for Matlab for me.
# IFLAGS, CFLAGS, OFLAGS  Flags to go with the compiler and linker.
#
# This Makefile works under Solaris 2.  To use it with a different Unix, 
# you may need to change the references to "ar" and/or "ranlib" and
# possibly other things.
#
# John Gilbert, 1996.
# Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
# See file ../copyright.m for complete copyright and licensing notice.

#MATLAB =        /import/matlab
MATLAB =        /usr/local/libexec/matlab

#CHACO =         /project/sparse/chaco/version2/code
CHACO =         ../../Chaco-2.0/code

MLCHACO =       .
DEST_DIR =      ..
#CC =            gcc 
CC =            cc 
IFLAGS =	-I${MATLAB}/extern/include -I${CHACO}/main
CFLAGS =        -Xa -G -xO4 -xcg92
OFLAGS =        -xO4 -xcg92
AR =             /usr/ccs/bin/ar rcv   # for solaris 2
RANLIB =         /usr/ccs/bin/ranlib   # for solaris 2


MLFILES.c=	${MLCHACO}/mlchaco.c \
		${MLCHACO}/bail.c \
		${MLCHACO}/chaco_check_graph.c \
		${MLCHACO}/check_input.c \
                ${MLCHACO}/smalloc.c \
                ${MLCHACO}/user_params.c

CHFILES.c=	${CHACO}/main/interface.c \
		${CHACO}/main/main.c \
		${CHACO}/submain/balance.c \
		${CHACO}/submain/divide.c \
		${CHACO}/submain/submain.c \
		${CHACO}/input/input_assign.c \
		${CHACO}/input/input.c \
		${CHACO}/input/input_geom.c \
		${CHACO}/input/input_graph.c \
		${CHACO}/input/read_params.c \
		${CHACO}/input/reflect_input.c \
		${CHACO}/input/read_val.c \
		${CHACO}/graph/free_graph.c \
		${CHACO}/graph/reformat.c \
		${CHACO}/graph/subgraph.c \
		${CHACO}/graph/graph_out.c \
		${CHACO}/inertial/eigenvec2.c \
		${CHACO}/inertial/eigenvec3.c \
		${CHACO}/inertial/inertial.c \
		${CHACO}/inertial/inertial1d.c \
		${CHACO}/inertial/inertial2d.c \
		${CHACO}/inertial/inertial3d.c \
		${CHACO}/inertial/make_subgeom.c \
		${CHACO}/klspiff/buckets.c \
		${CHACO}/klspiff/buckets_bi.c \
		${CHACO}/klspiff/buckets1.c \
		${CHACO}/klspiff/bilistops.c \
		${CHACO}/klspiff/coarsen_kl.c \
		${CHACO}/klspiff/count_weights.c \
		${CHACO}/klspiff/compress_ewgts.c \
		${CHACO}/klspiff/kl_init.c \
		${CHACO}/klspiff/kl_output.c \
		${CHACO}/klspiff/klspiff.c \
		${CHACO}/klspiff/make_bndy_list.c \
		${CHACO}/klspiff/make_kl_list.c \
		${CHACO}/klspiff/nway_kl.c \
		${CHACO}/klvspiff/bpm_improve.c \
		${CHACO}/klvspiff/bucketsv.c \
		${CHACO}/klvspiff/clear_dvals.c \
		${CHACO}/klvspiff/coarsen_klv.c \
		${CHACO}/klvspiff/countup_vtx_sep.c \
		${CHACO}/klvspiff/find_bndy.c \
		${CHACO}/klvspiff/klv_init.c \
		${CHACO}/klvspiff/klvspiff.c \
		${CHACO}/klvspiff/make_bpgraph.c \
		${CHACO}/klvspiff/make_sep_list.c \
		${CHACO}/klvspiff/matching.c \
		${CHACO}/klvspiff/nway_klv.c \
		${CHACO}/coarsen/coarsen.c \
		${CHACO}/coarsen/interpolate.c \
		${CHACO}/coarsen/makecgraph.c \
		${CHACO}/coarsen/coarsen1.c \
		${CHACO}/coarsen/countcedges.c \
		${CHACO}/coarsen/makev2cv.c \
		${CHACO}/coarsen/maxmatch.c \
		${CHACO}/coarsen/maxmatch1.c \
		${CHACO}/coarsen/maxmatch2.c  \
		${CHACO}/coarsen/maxmatch3.c \
		${CHACO}/coarsen/maxmatch4.c \
		${CHACO}/connect/add_edges.c \
		${CHACO}/connect/connected.c \
		${CHACO}/connect/find_edges.c \
		${CHACO}/eigen/bidir.c \
		${CHACO}/eigen/bisect.c \
		${CHACO}/eigen/checkeig.c \
		${CHACO}/eigen/checkeig_ext.c \
		${CHACO}/eigen/checkorth.c \
		${CHACO}/eigen/cksturmcnt.c \
		${CHACO}/eigen/mkeigvecs.c\
		${CHACO}/eigen/eigensolve.c \
		${CHACO}/eigen/get_extval.c \
		${CHACO}/eigen/get_ritzvals.c \
		${CHACO}/eigen/lanczos_FO.c \
		${CHACO}/eigen/lanczos_SO.c \
		${CHACO}/eigen/lanczos_SO_float.c \
		${CHACO}/eigen/lanczos_ext.c \
		${CHACO}/eigen/lanczos_ext_float.c \
		${CHACO}/eigen/lanc_seconds.c\
		${CHACO}/eigen/lanpause.c \
		${CHACO}/eigen/makeorthlnk.c \
		${CHACO}/eigen/mkscanlist.c \
		${CHACO}/eigen/orthog1.c \
		${CHACO}/eigen/orthogonalize.c \
		${CHACO}/eigen/orthogvec.c \
		${CHACO}/eigen/ql.c \
		${CHACO}/eigen/rqi.c \
		${CHACO}/eigen/rqi_ext.c \
		${CHACO}/eigen/scale_diag.c \
		${CHACO}/eigen/scanmax.c \
		${CHACO}/eigen/scanmin.c \
		${CHACO}/eigen/solistout.c \
		${CHACO}/eigen/sorthog.c \
		${CHACO}/eigen/splarax.c \
		${CHACO}/eigen/sturmcnt.c \
		${CHACO}/eigen/Tevec.c \
		${CHACO}/eigen/tri_solve.c \
		${CHACO}/eigen/warnings.c \
		${CHACO}/symmlq/aprod.c \
		${CHACO}/symmlq/msolve.c \
		${CHACO}/symmlq/pow_dd.c \
		${CHACO}/symmlq/symmlq.c \
		${CHACO}/symmlq/symmlqblas.c \
		${CHACO}/tinvit/tinvit.c \
		${CHACO}/tinvit/pythag.c \
		${CHACO}/tinvit/epslon.c \
		${CHACO}/optimize/determinant.c \
		${CHACO}/optimize/func2d.c \
		${CHACO}/optimize/func3d.c \
		${CHACO}/optimize/opt2d.c \
		${CHACO}/optimize/opt3d.c \
		${CHACO}/assign/assign.c \
		${CHACO}/assign/assign_out.c \
		${CHACO}/assign/mapper.c \
		${CHACO}/assign/median.c \
		${CHACO}/assign/merge_assign.c \
		${CHACO}/assign/rec_median.c \
		${CHACO}/assign/rotate.c \
		${CHACO}/assign/y2x.c \
		${CHACO}/bpmatch/checkbp.c \
		${CHACO}/bpmatch/inits2d.c \
		${CHACO}/bpmatch/inits3d.c \
		${CHACO}/bpmatch/genvals2d.c \
		${CHACO}/bpmatch/genvals3d.c \
		${CHACO}/bpmatch/map2d.c \
		${CHACO}/bpmatch/map3d.c \
		${CHACO}/bpmatch/movevtxs.c \
		${CHACO}/bpmatch/sorts2d.c \
		${CHACO}/bpmatch/sorts3d.c \
		${CHACO}/refine_map/compute_cube_edata.c \
		${CHACO}/refine_map/compute_cube_vdata.c \
		${CHACO}/refine_map/refine_cube.c \
		${CHACO}/refine_map/update_cube_edata.c \
		${CHACO}/refine_map/update_cube_vdata.c \
		${CHACO}/refine_map/find_edge_cube.c \
		${CHACO}/refine_map/init_cube_edata.c \
		${CHACO}/refine_map/compute_mesh_edata.c \
		${CHACO}/refine_map/compute_mesh_vdata.c \
		${CHACO}/refine_map/find_edge_mesh.c \
		${CHACO}/refine_map/init_mesh_edata.c \
		${CHACO}/refine_map/refine_mesh.c \
		${CHACO}/refine_map/update_mesh_edata.c \
		${CHACO}/refine_map/update_mesh_vdata.c \
		${CHACO}/refine_map/refine_map.c \
		${CHACO}/refine_map/make_comm_graph.c \
		${CHACO}/refine_part/refine_part.c \
		${CHACO}/refine_part/kl_refine.c \
		${CHACO}/refine_part/make_maps_ref.c \
		${CHACO}/refine_part/make_terms_ref.c \
		${CHACO}/internal/force_internal.c \
		${CHACO}/internal/improve_internal.c \
		${CHACO}/internal/check_internal.c \
		${CHACO}/misc/define_subcubes.c \
		${CHACO}/misc/define_submeshes.c \
		${CHACO}/misc/divide_procs.c \
		${CHACO}/misc/merge_goals.c \
		${CHACO}/misc/make_term_props.c \
		${CHACO}/misc/count.c \
		${CHACO}/misc/countup.c \
		${CHACO}/misc/countup_cube.c \
		${CHACO}/misc/countup_mesh.c \
		${CHACO}/misc/make_subgoal.c \
		${CHACO}/misc/find_maxdeg.c \
		${CHACO}/misc/make_maps.c \
		${CHACO}/misc/make_setlists.c \
		${CHACO}/misc/sequence.c \
		${CHACO}/misc/perturb.c \
		${CHACO}/misc/simple_part.c \
		${CHACO}/misc/time_kernels.c \
		${CHACO}/misc/timing.c \
		${CHACO}/util/affirm.c \
		${CHACO}/util/array_alloc_2D.c \
		${CHACO}/util/bit_reverse.c \
		${CHACO}/util/checkpnt.c \
		${CHACO}/util/cpvec.c \
		${CHACO}/util/dot.c \
		${CHACO}/util/doubleout.c \
		${CHACO}/util/input_int.c \
		${CHACO}/util/gray.c \
		${CHACO}/util/machine_params.c \
		${CHACO}/util/makevwsqrt.c \
		${CHACO}/util/mkvec.c \
		${CHACO}/util/norm.c \
		${CHACO}/util/normalize.c \
		${CHACO}/util/mergesort.c \
		${CHACO}/util/randomize.c \
		${CHACO}/util/scadd.c \
		${CHACO}/util/seconds.c \
		${CHACO}/util/setvec.c \
		${CHACO}/util/shell_sort.c \
		${CHACO}/util/strout.c \
		${CHACO}/util/tri_prod.c \
		${CHACO}/util/true_or_false.c \
		${CHACO}/util/update.c  \
		${CHACO}/util/vecout.c \
		${CHACO}/util/vecran.c \
		${CHACO}/util/vecscale.c 

MLFILES.o=      $(MLFILES.c:.c=.o)

CHFILES.o=      $(CHFILES.c:.c=.o)

#mlchaco:        ${MLFILES.o} chaco.a Makefile
#		${MATLAB}/bin/cmex CC='gcc -G' -lm ${OFLAGS} ${MLFILES.o} chaco.a; \
#                 mv mlchaco.mex* ${DEST_DIR}

mlchaco:	${MLFILES.c} chaco.a Makefile
		mex -V4 -output mlchaco ${MLFILES.c} chaco.a -I${CHACO}/main
		mv mlchaco.mex* ${DEST_DIR}

chaco.a:        ${CHFILES.o}
		${AR} chaco.a ${CHFILES.o} ; ${RANLIB} chaco.a

clean:
		rm -f ${CHACO}/*/*.o chaco.a ${MLCHACO}/*.o ${DEST_DIR}/core

.c.o:
		${CC} -c ${IFLAGS} ${CFLAGS} -o $*.o $*.c

