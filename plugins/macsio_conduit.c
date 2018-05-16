/*
Copyright (c) 2018, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Brad Whitlock

LLNL-CODE-676051. All rights reserved.

This file is part of MACSio

Please also read the LICENSE file at the top of the source code directory or
folder hierarchy.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License (as published by the Free Software
Foundation) version 2, dated June 1991.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <json-cwx/json.h>

#include <macsio_clargs.h>
#include <macsio_iface.h>
#include <macsio_log.h>
#include <macsio_main.h>
#include <macsio_mif.h>
#include <macsio_utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <conduit/conduit.h>
#include <conduit/conduit_blueprint.h>
#include <conduit/conduit_relay.h>
#include <conduit/conduit_utils.h>

#define MAXKEYLEN 200

/* Disable debugging messages */

/*!
\addtogroup plugins
@{
*/

/*!
\addtogroup Conduit
@{
*/

/* the name you want to assign to the interface */
static char const *iface_name = "conduit";
static char const *iface_ext = "blueprint_root";

#if 0 /* HDF5 Junk */
static int use_log = 0;
static int no_collective = 0;
static int no_single_chunk = 0;
static int silo_block_size = 0;
static int silo_block_count = 0;
static int sbuf_size = -1;
static int mbuf_size = -1;
static int rbuf_size = -1;
static int lbuf_size = 0;
static const char *filename;
static hid_t fid;
static hid_t dspc = -1;
static int show_errors = 0;
static char compression_alg_str[64];
static char compression_params_str[512];

static hid_t make_fapl()
{
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    herr_t h5status = 0;

    if (sbuf_size >= 0)
        h5status |= H5Pset_sieve_buf_size(fapl_id, sbuf_size);

    if (mbuf_size >= 0)
        h5status |= H5Pset_meta_block_size(fapl_id, mbuf_size);

    if (rbuf_size >= 0)
        h5status |= H5Pset_small_data_block_size(fapl_id, mbuf_size);

#if 0
    if (silo_block_size && silo_block_count)
    {
        h5status |= H5Pset_fapl_silo(fapl_id);
        h5status |= H5Pset_silo_block_size_and_count(fapl_id, (hsize_t) silo_block_size,
            silo_block_count);
    }
    else
    if (use_log)
    {
        int flags = H5FD_LOG_LOC_IO|H5FD_LOG_NUM_IO|H5FD_LOG_TIME_IO|H5FD_LOG_ALLOC;

        if (lbuf_size > 0)
            flags = H5FD_LOG_ALL;

        h5status |= H5Pset_fapl_log(fapl_id, "macsio_hdf5_log.out", flags, lbuf_size);
    }
#endif

    {
        H5AC_cache_config_t config;

        /* Acquire a default mdc config struct */
        config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
        H5Pget_mdc_config(fapl_id, &config);
#define MAINZER_PARAMS 1
#if MAINZER_PARAMS
        config.set_initial_size = (hbool_t) 1;
        config.initial_size = 16 * 1024;
        config.min_size = 8 * 1024;
        config.epoch_length = 3000;
        config.lower_hr_threshold = 0.95;
#endif
        H5Pset_mdc_config(fapl_id, &config);
    }

    if (h5status < 0)
    {
        if (fapl_id >= 0)
            H5Pclose(fapl_id);
        return 0;
    }

    return fapl_id;
}

static int
get_tokval(char const *src_str, char const *token_to_match, void *val_ptr)
{
    int toklen;
    char dummy[16];
    void *val_ptr_tmp = &dummy[0];

    if (!src_str) return 0;
    if (!token_to_match) return 0;

    toklen = strlen(token_to_match)-2;

    if (strncasecmp(src_str, token_to_match, toklen))
        return 0;

    /* First, check matching sscanf *without* risk of writing to val_ptr */
    if (sscanf(&src_str[toklen], &token_to_match[toklen], val_ptr_tmp) != 1)
        return 0;

    sscanf(&src_str[toklen], &token_to_match[toklen], val_ptr);
    return 1;
}

static hid_t make_dcpl(char const *alg_str, char const *params_str, hid_t space_id, hid_t dtype_id)
{
    int shuffle = -1;
    int minsize = -1;
    int gzip_level = -1;
    int zfp_precision = -1;
    unsigned int szip_pixels_per_block = 0;
    float zfp_rate = -1;
    float zfp_accuracy = -1;
    char szip_method[64], szip_chunk_str[64];
    char *token, *string, *tofree;
    int ndims;
    hsize_t dims[4], maxdims[4];
    hid_t retval = H5Pcreate(H5P_DATASET_CREATE);

    szip_method[0] = '\0';
    szip_chunk_str[0] = '\0';

    /* Initially, set contiguous layout. May reset to chunked later */
    H5Pset_layout(retval, H5D_CONTIGUOUS);

    if (!alg_str || !strlen(alg_str))
        return retval;

    ndims = H5Sget_simple_extent_ndims(space_id);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);

    /* We can make a pass through params string without being specific about
       algorithm because there are presently no symbol collisions there */
    tofree = string = strdup(params_str);
    while ((token = strsep(&string, ",")) != NULL)
    {
        if (get_tokval(token, "minsize=%d", &minsize))
            continue;
        if (get_tokval(token, "shuffle=%d", &shuffle))
            continue;
        if (get_tokval(token, "level=%d", &gzip_level))
            continue;
        if (get_tokval(token, "rate=%f", &zfp_rate))
            continue;
        if (get_tokval(token, "precision=%d", &zfp_precision))
            continue;
        if (get_tokval(token, "accuracy=%f", &zfp_accuracy))
            continue;
        if (get_tokval(token, "method=%s", szip_method))
            continue;
        if (get_tokval(token, "block=%u", &szip_pixels_per_block))
            continue;
        if (get_tokval(token, "chunk=%s", szip_chunk_str))
            continue;
    }
    free(tofree);

    /* check for minsize compression threshold */
    minsize = minsize != -1 ? minsize : 1024;
    if (H5Sget_simple_extent_npoints(space_id) < minsize)
        return retval;

    /*
     * Ok, now handle various properties related to compression
     */
 
    /* Initially, as a default in case nothing else is selected,
       set chunk size equal to dataset size (e.g. single chunk) */
    H5Pset_chunk(retval, ndims, dims);

    if (!strncasecmp(alg_str, "gzip", 4))
    {
        if (shuffle == -1 || shuffle == 1)
            H5Pset_shuffle(retval);
        H5Pset_deflate(retval, gzip_level!=-1?gzip_level:9);
    }
    else if (!strncasecmp(alg_str, "szip", 4))
    {
        static int have_issued_warning = 0;
        if (!have_issued_warning)
            MACSIO_LOG_MSG(Warn, ("szip compressor not available in this build"));
        have_issued_warning = 1;
    }

    return retval;
}

static void main_dump_sif(json_object *main_obj, int dumpn, double dumpt)
{
#ifdef HAVE_MPI
    int ndims;
    int i, v, p;
    char const *mesh_type = json_object_path_get_string(main_obj, "clargs/part_type");
    char fileName[256];
    int use_part_count;

    hid_t h5file_id;
    hid_t fapl_id = make_fapl();
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t null_space_id = H5Screate(H5S_NULL);
    hid_t fspace_nodal_id, fspace_zonal_id;
    hsize_t global_log_dims_nodal[3];
    hsize_t global_log_dims_zonal[3];

    MPI_Info mpiInfo = MPI_INFO_NULL;

//#warning WE ARE DOING SIF SLIGHTLY WRONG, DUPLICATING SHARED NODES
//#warning INCLUDE ARGS FOR ISTORE AND K_SYM
//#warning INCLUDE ARG PROCESS FOR HINTS
//#warning FAPL PROPS: ALIGNMENT 
#if H5_HAVE_PARALLEL
    H5Pset_fapl_mpio(fapl_id, MACSIO_MAIN_Comm, mpiInfo);
#endif

//#warning FOR MIF, NEED A FILEROOT ARGUMENT OR CHANGE TO FILEFMT ARGUMENT
    /* Construct name for the HDF5 file */
    sprintf(fileName, "%s_hdf5_%03d.%s",
        json_object_path_get_string(main_obj, "clargs/filebase"),
        dumpn,
        json_object_path_get_string(main_obj, "clargs/fileext"));

    MACSIO_UTILS_RecordOutputFiles(dumpn, fileName);
    h5file_id = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

    /* Create an HDF5 Dataspace for the global whole of mesh and var objects in the file. */
    ndims = json_object_path_get_int(main_obj, "clargs/part_dim");
    json_object *global_log_dims_array =
        json_object_path_get_array(main_obj, "problem/global/LogDims");
    json_object *global_parts_log_dims_array =
        json_object_path_get_array(main_obj, "problem/global/PartsLogDims");
    /* Note that global zonal array is smaller in each dimension by one *ON*EACH*BLOCK*
       in the associated dimension. */
    for (i = 0; i < ndims; i++)
    {
        int parts_log_dims_val = JsonGetInt(global_parts_log_dims_array, "", i);
        global_log_dims_nodal[ndims-1-i] = (hsize_t) JsonGetInt(global_log_dims_array, "", i);
        global_log_dims_zonal[ndims-1-i] = global_log_dims_nodal[ndims-1-i] -
            JsonGetInt(global_parts_log_dims_array, "", i);
    }
    fspace_nodal_id = H5Screate_simple(ndims, global_log_dims_nodal, 0);
    fspace_zonal_id = H5Screate_simple(ndims, global_log_dims_zonal, 0);

    /* Get the list of vars on the first part as a guide to loop over vars */
    json_object *part_array = json_object_path_get_array(main_obj, "problem/parts");
    json_object *first_part_obj = json_object_array_get_idx(part_array, 0);
    json_object *first_part_vars_array = json_object_path_get_array(first_part_obj, "Vars");

    /* Dataset transfer property list used in all H5Dwrite calls */
#if H5_HAVE_PARALLEL
    if (no_collective)
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
    else
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
#endif


    /* Loop over vars and then over parts */
    /* currently assumes all vars exist on all ranks. but not all parts */
    for (v = -1; v < json_object_array_length(first_part_vars_array); v++) /* -1 start is for Mesh */
    {

//#warning SKIPPING MESH
        if (v == -1) continue; /* All ranks skip mesh (coords) for now */

        /* Inspect the first part's var object for name, datatype, etc. */
        json_object *var_obj = json_object_array_get_idx(first_part_vars_array, v);
        char const *varName = json_object_path_get_string(var_obj, "name");
        char *centering = strdup(json_object_path_get_string(var_obj, "centering"));
        json_object *dataobj = json_object_path_get_extarr(var_obj, "data");
//#warning JUST ASSUMING TWO TYPES NOW. CHANGE TO A FUNCTION
        hid_t dtype_id = json_object_extarr_type(dataobj)==json_extarr_type_flt64? 
                H5T_NATIVE_DOUBLE:H5T_NATIVE_INT;
        hid_t fspace_id = H5Scopy(strcmp(centering, "zone") ? fspace_nodal_id : fspace_zonal_id);
        hid_t dcpl_id = make_dcpl(compression_alg_str, compression_params_str, fspace_id, dtype_id);

        /* Create the file dataset (using old-style H5Dcreate API here) */
//#warning USING DEFAULT DCPL: LATER ADD COMPRESSION, ETC.
        
        hid_t ds_id = H5Dcreate1(h5file_id, varName, dtype_id, fspace_id, dcpl_id); 
        H5Sclose(fspace_id);
        H5Pclose(dcpl_id);

        /* Loop to make write calls for this var for each part on this rank */
//#warning USE NEW MULTI-DATASET API WHEN AVAILABLE TO AGLOMERATE ALL PARTS INTO ONE CALL
        use_part_count = (int) ceil(json_object_path_get_double(main_obj, "clargs/avg_num_parts"));
        for (p = 0; p < use_part_count; p++)
        {
            json_object *part_obj = json_object_array_get_idx(part_array, p);
            json_object *var_obj = 0;
            hid_t mspace_id = H5Scopy(null_space_id);
            void const *buf = 0;

            fspace_id = H5Scopy(null_space_id);

            /* this rank actually has something to contribute to the H5Dwrite call */
            if (part_obj)
            {
                int i;
                hsize_t starts[3], counts[3];
                json_object *vars_array = json_object_path_get_array(part_obj, "Vars");
                json_object *mesh_obj = json_object_path_get_object(part_obj, "Mesh");
                json_object *var_obj = json_object_array_get_idx(vars_array, v);
                json_object *extarr_obj = json_object_path_get_extarr(var_obj, "data");
                json_object *global_log_origin_array =
                    json_object_path_get_array(part_obj, "GlobalLogOrigin");
                json_object *global_log_indices_array =
                    json_object_path_get_array(part_obj, "GlobalLogIndices");
                json_object *mesh_dims_array = json_object_path_get_array(mesh_obj, "LogDims");
                for (i = 0; i < ndims; i++)
                {
                    starts[ndims-1-i] =
                        json_object_get_int(json_object_array_get_idx(global_log_origin_array,i));
                    counts[ndims-1-i] =
                        json_object_get_int(json_object_array_get_idx(mesh_dims_array,i));
                    if (!strcmp(centering, "zone"))
                    {
                        counts[ndims-1-i]--;
                        starts[ndims-1-i] -=
                            json_object_get_int(json_object_array_get_idx(global_log_indices_array,i));
                    }
                }

                /* set selection of filespace */
                fspace_id = H5Dget_space(ds_id);
                H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, starts, 0, counts, 0);

                /* set dataspace of data in memory */
                mspace_id = H5Screate_simple(ndims, counts, 0);
                buf = json_object_extarr_data(extarr_obj);
            }

            H5Dwrite(ds_id, dtype_id, mspace_id, fspace_id, dxpl_id, buf);
            H5Sclose(fspace_id);
            H5Sclose(mspace_id);

        }

        H5Dclose(ds_id);
        free(centering);
    }

    H5Sclose(fspace_nodal_id);
    H5Sclose(fspace_zonal_id);
    H5Sclose(null_space_id);
    H5Pclose(dxpl_id);
    H5Pclose(fapl_id);
    H5Fclose(h5file_id);

#endif
}

typedef struct _user_data {
    hid_t groupId;
} user_data_t;

static void *CreateHDF5File(const char *fname, const char *nsname, void *userData)
{
    hid_t *retval = 0;
    hid_t h5File = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h5File >= 0)
    {
//#warning USE NEWER GROUP CREATION SETTINGS OF HDF5
        if (nsname && userData)
        {
            user_data_t *ud = (user_data_t *) userData;
            ud->groupId = H5Gcreate1(h5File, nsname, 0);
        }
        retval = (hid_t *) malloc(sizeof(hid_t));
        *retval = h5File;
    }
    return (void *) retval;
}

static void *OpenHDF5File(const char *fname, const char *nsname,
                   MACSIO_MIF_ioFlags_t ioFlags, void *userData)
{
    hid_t *retval;
    hid_t h5File = H5Fopen(fname, ioFlags.do_wr ? H5F_ACC_RDWR : H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5File >= 0)
    {
        if (ioFlags.do_wr && nsname && userData)
        {
            user_data_t *ud = (user_data_t *) userData;
            ud->groupId = H5Gcreate1(h5File, nsname, 0);
        }
        retval = (hid_t *) malloc(sizeof(hid_t));
        *retval = h5File;
    }
    return (void *) retval;
}

static void CloseHDF5File(void *file, void *userData)
{
    if (userData)
    {
        user_data_t *ud = (user_data_t *) userData;
        H5Gclose(ud->groupId);
    }
    H5Fclose(*((hid_t*) file));
    free(file);
}

static void write_mesh_part(hid_t h5loc, json_object *part_obj)
{
//#warning WERE SKPPING THE MESH (COORDS) OBJECT PRESENTLY
    int i;
    json_object *vars_array = json_object_path_get_array(part_obj, "Vars");

    for (i = 0; i < json_object_array_length(vars_array); i++)
    {
        int j;
        hsize_t var_dims[3];
        hid_t fspace_id, ds_id, dcpl_id;
        json_object *var_obj = json_object_array_get_idx(vars_array, i);
        json_object *data_obj = json_object_path_get_extarr(var_obj, "data");
        char const *varname = json_object_path_get_string(var_obj, "name");
        int ndims = json_object_extarr_ndims(data_obj);
        void const *buf = json_object_extarr_data(data_obj);
        hid_t dtype_id = json_object_extarr_type(data_obj)==json_extarr_type_flt64? 
                H5T_NATIVE_DOUBLE:H5T_NATIVE_INT;

        for (j = 0; j < ndims; j++)
            var_dims[j] = json_object_extarr_dim(data_obj, j);

        fspace_id = H5Screate_simple(ndims, var_dims, 0);
        dcpl_id = make_dcpl(compression_alg_str, compression_params_str, fspace_id, dtype_id);
        ds_id = H5Dcreate1(h5loc, varname, dtype_id, fspace_id, dcpl_id); 
        H5Dwrite(ds_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(ds_id);
        H5Pclose(dcpl_id);
        H5Sclose(fspace_id);
    }
}

static void main_dump_mif(json_object *main_obj, int numFiles, int dumpn, double dumpt)
{
    int size, rank;
    hid_t *h5File_ptr;
    hid_t h5File;
    hid_t h5Group;
    char fileName[256];
    int i, len;
    int *theData;
    user_data_t userData;
    MACSIO_MIF_ioFlags_t ioFlags = {MACSIO_MIF_WRITE,
        JsonGetInt(main_obj, "clargs/exercise_scr")&0x1};

//#warning MAKE WHOLE FILE USE HDF5 1.8 INTERFACE
//#warning SET FILE AND DATASET PROPERTIES
//#warning DIFFERENT MPI TAGS FOR DIFFERENT PLUGINS AND CONTEXTS
    MACSIO_MIF_baton_t *bat = MACSIO_MIF_Init(numFiles, ioFlags, MACSIO_MAIN_Comm, 3,
        CreateHDF5File, OpenHDF5File, CloseHDF5File, &userData);

    rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
    size = json_object_path_get_int(main_obj, "parallel/mpi_size");

    /* Construct name for the silo file */
    sprintf(fileName, "%s_hdf5_%05d_%03d.%s",
        json_object_path_get_string(main_obj, "clargs/filebase"),
        MACSIO_MIF_RankOfGroup(bat, rank),
        dumpn,
        json_object_path_get_string(main_obj, "clargs/fileext"));

    MACSIO_UTILS_RecordOutputFiles(dumpn, fileName);
    
    h5File_ptr = (hid_t *) MACSIO_MIF_WaitForBaton(bat, fileName, 0);
    h5File = *h5File_ptr;
    h5Group = userData.groupId;

    json_object *parts = json_object_path_get_array(main_obj, "problem/parts");

    for (int i = 0; i < json_object_array_length(parts); i++)
    {
        char domain_dir[256];
        json_object *this_part = json_object_array_get_idx(parts, i);
        hid_t domain_group_id;

        snprintf(domain_dir, sizeof(domain_dir), "domain_%07d",
            json_object_path_get_int(this_part, "Mesh/ChunkID"));
 
        domain_group_id = H5Gcreate1(h5File, domain_dir, 0);

        write_mesh_part(domain_group_id, this_part);

        H5Gclose(domain_group_id);
    }

    /* If this is the 'root' processor, also write Silo's multi-XXX objects */
#if 0
    if (rank == 0)
        WriteMultiXXXObjects(main_obj, siloFile, bat);
#endif

    /* Hand off the baton to the next processor. This winds up closing
     * the file so that the next processor that opens it can be assured
     * of getting a consistent and up to date view of the file's contents. */
    MACSIO_MIF_HandOffBaton(bat, h5File_ptr);

    /* We're done using MACSIO_MIF, so finish it off */
    MACSIO_MIF_Finish(bat);

}

static void main_dump(int argi, int argc, char **argv, json_object *main_obj,
    int dumpn, double dumpt)
{
    int rank, size, numFiles;

//#warning SET ERROR MODE OF HDF5 LIBRARY

    /* Without this barrier, I get strange behavior with Silo's MACSIO_MIF interface */
    mpi_errno = MPI_Barrier(MACSIO_MAIN_Comm);

    /* process cl args */
    process_args(argi, argc, argv);

    rank = json_object_path_get_int(main_obj, "parallel/mpi_rank");
    size = json_object_path_get_int(main_obj, "parallel/mpi_size");

//#warning MOVE TO A FUNCTION
    /* ensure we're in MIF mode and determine the file count */
    json_object *parfmode_obj = json_object_path_get_array(main_obj, "clargs/parallel_file_mode");
    if (parfmode_obj)
    {
        json_object *modestr = json_object_array_get_idx(parfmode_obj, 0);
        json_object *filecnt = json_object_array_get_idx(parfmode_obj, 1);
//#warning ERRORS NEED TO GO TO LOG FILES AND ERROR BEHAVIOR NEEDS TO BE HONORED
        if (!strcmp(json_object_get_string(modestr), "SIF"))
        {
            main_dump_sif(main_obj, dumpn, dumpt);
        }
        else
        {
            numFiles = json_object_get_int(filecnt);
            main_dump_mif(main_obj, numFiles, dumpn, dumpt);
        }
    }
    else
    {
        char const * modestr = json_object_path_get_string(main_obj, "clargs/parallel_file_mode");
        if (!strcmp(modestr, "SIF"))
        {
            float avg_num_parts = json_object_path_get_double(main_obj, "clargs/avg_num_parts");
            if (avg_num_parts == (float ((int) avg_num_parts)))
                main_dump_sif(main_obj, dumpn, dumpt);
            else
            {
//#warning CURRENTLY, SIF CAN WORK ONLY ON WHOLE PART COUNTS
                MACSIO_LOG_MSG(Die, ("HDF5 plugin cannot currently handle SIF mode where "
                    "there are different numbers of parts on each MPI rank. "
                    "Set --avg_num_parts to an integral value." ));
            }
        }
        else if (!strcmp(modestr, "MIFMAX"))
            numFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
        else if (!strcmp(modestr, "MIFAUTO"))
        {
            /* Call utility to determine optimal file count */
//#warning ADD UTILIT TO DETERMINE OPTIMAL FILE COUNT
        }
        main_dump_mif(main_obj, numFiles, dumpn, dumpt);
    }
}
#endif /*HDF5 Junk*/

/* Static vars */
static conduit_node *g_conduit_about = NULL;
static conduit_node *g_conduit_relay_about = NULL;
static char g_preferred_protocol[100];
static char g_preferred_protocol_ext[100];

/*-----------------------------------------------------------------------------*/
static void
set_preferred_protocol(void)
{
#if 0
    strcpy(g_preferred_protocol, "json");
    strcpy(g_preferred_protocol_ext, ".json");
    return;
#endif

#if 0
    if(conduit_node_has_path(g_conduit_relay_about, "io/protocols/conduit_silo_mesh"))
    {
        strcpy(g_preferred_protocol, "conduit_silo_mesh");
        strcpy(g_preferred_protocol_ext, ".silo");
        return;
    }
#endif

    if(conduit_node_has_path(g_conduit_relay_about, "io/protocols/hdf5"))
    {
        strcpy(g_preferred_protocol, "hdf5");
        strcpy(g_preferred_protocol_ext, ".h5");
        return;
    }

    if(conduit_node_has_path(g_conduit_relay_about, "io/protocols/adios"))
    {
        strcpy(g_preferred_protocol, "adios");
        strcpy(g_preferred_protocol_ext, ".bp");
        /* Let the user pick a transport too (bp, hdf5, netcdf, flexpath, ...) */
        return;
    }

    strcpy(g_preferred_protocol, "conduit_bin");
    strcpy(g_preferred_protocol_ext, ".bin");
}

/*-----------------------------------------------------------------------------*/
/* This is an entry point from MACSio to the conduit plugin. */
static int process_args(int argi, int argc, char *argv[])
{
    const MACSIO_CLARGS_ArgvFlags_t argFlags = {MACSIO_CLARGS_WARN, MACSIO_CLARGS_TOMEM};

    char *c_protocol = NULL;

#if 1
    /* This thing seems to let MACSio parse args and put values into static global vars in this driver. */
    MACSIO_CLARGS_ProcessCmdline(0, argFlags, argi, argc, argv,
        "--protocol %s", "conduit_silo_mesh", /*MACSIO_CLARGS_NODEFAULT,*/
            "Output protocol for Conduit",
            &c_protocol,
           MACSIO_CLARGS_END_OF_ARGS);
#else
    MACSIO_CLARGS_ProcessCmdline(0, argFlags, argi, argc, argv,
        "--show_errors", "",
            "Show low-level HDF5 errors",
            &c_format,
        "--compression %s %s", MACSIO_CLARGS_NODEFAULT,
            "The first string argument is the compression algorithm name. The second\n"
            "string argument is a comma-separated set of params of the form\n"
            "'param1=val1,param2=val2,param3=val3. The various algorithm names and\n"
            "their parameter meanings are described below. Note that some parameters are\n"
            "not specific to any algorithm. Those are described first followed by\n"
            "individual algorithm-specific parameters.\n"
            "\n"
            "minsize=%d : min. size of dataset (in terms of a count of values)\n"
            "    upon which compression will even be attempted. Default is 1024.\n"
            "shuffle=<int>: Boolean (zero or non-zero) to indicate whether to use\n"
            "    HDF5's byte shuffling filter *prior* to compression. Default depends\n"
            "    on algorithm. By default, shuffling is NOT used for lindstrom-zfp but IS\n"
            "    used with all other algorithms.\n"
            "\n"

            "\"szip\"\n"
            "    method=%s : specify 'ec' for entropy coding or 'nn' for nearest\n"
            "        neighbor. Default is 'nn'\n"
            "    block=%d : (pixels-per-block) must be an even integer <= 32. See\n"
            "        See H5Pset_szip in HDF5 documentation for more information.\n"
            "        Default is 32.\n"
            "    chunk=%d:%d : colon-separated dimensions specifying chunk size in\n"
            "        each dimension higher than the first (fastest varying) dimension.\n"
            "\n"
            "\"gzip\"\n"
            "    level=%d : A value in the range [1,9], inclusive, trading off time to\n"
            "        compress with amount of compression. Level=1 results in best speed\n"
            "        but worst compression whereas level=9 results in best compression\n"
            "        but worst speed. Values outside [1,9] are clamped. Default is 9.\n"
            "\n"
            "Examples:\n"
            "    --compression lindstrom-zfp rate=18.5\n"
            "    --compression gzip minsize=1024,level=9\n"
            "    --compression szip shuffle=0,options=nn,pixels_per_block=16\n"
            "\n",
            &c_alg, &c_params,
        "--no_collective", "",
            "Use independent, not collective, I/O calls in SIF mode.",
            &no_collective,
        "--no_single_chunk", "",
            "Do not single chunk the datasets (currently ignored).",
            &no_single_chunk,
        "--sieve_buf_size %d", MACSIO_CLARGS_NODEFAULT,
            "Specify sieve buffer size (see H5Pset_sieve_buf_size)",
            &sbuf_size,
        "--meta_block_size %d", MACSIO_CLARGS_NODEFAULT,
            "Specify size of meta data blocks (see H5Pset_meta_block_size)",
            &mbuf_size,
        "--small_block_size %d", MACSIO_CLARGS_NODEFAULT,
            "Specify threshold size for data blocks considered to be 'small'\n"
            "(see H5Pset_small_data_block_size)",
            &rbuf_size,
        "--log", "",
            "Use logging Virtual File Driver (see H5Pset_fapl_log)",
            &use_log,
           MACSIO_CLARGS_END_OF_ARGS);
#endif

    if(c_protocol != NULL)
        strcpy(g_preferred_protocol, c_protocol);

    /* Get a file format extension for the protocol -- this really ought to come from Conduit. */
    if(strcmp(g_preferred_protocol, "conduit_bin") == 0)
        strcpy(g_preferred_protocol_ext, ".bin");
    else if(strcmp(g_preferred_protocol, "json") == 0)
        strcpy(g_preferred_protocol_ext, ".json");
    else if(strcmp(g_preferred_protocol, "hdf5") == 0)
        strcpy(g_preferred_protocol_ext, ".h5");
    else if(strcmp(g_preferred_protocol, "conduit_silo") == 0)
        strcpy(g_preferred_protocol_ext, ".silo");
    else if(strcmp(g_preferred_protocol, "conduit_silo_mesh") == 0)
        strcpy(g_preferred_protocol_ext, ".silo");
    else if(strcmp(g_preferred_protocol, "mpi") == 0)
        strcpy(g_preferred_protocol_ext, ".mpi");
    else if(strcmp(g_preferred_protocol, "adios") == 0) /* Do we want to do "adios_<transport>"? */
        strcpy(g_preferred_protocol_ext, ".bp");
    else
    {
        set_preferred_protocol();
    }

    return 0;
}

/*-----------------------------------------------------------------------------*/
static void
json_object_to_blueprint_uniform(json_object *part, conduit_node *mesh,
    const char *meshName, const char *topoName)
{
    char key[MAXKEYLEN];

    /* Create the coordinates for a rectilinear mesh, zero copy. */
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "uniform");

    if(JsonGetInt(part, "Mesh/Coords/NumX") > 0)
    {
        conduit_node_set_path_int(mesh, "coordsets/coords/dims/i", JsonGetInt(part, "Mesh/Coords/NumX"));
        conduit_node_set_path_double(mesh, "coordsets/coords/origin/x", JsonGetDbl(part, "Mesh/Coords/OriginX"));
        conduit_node_set_path_double(mesh, "coordsets/coords/spacing/dx", JsonGetDbl(part, "Mesh/Coords/DeltaX"));
    }

    if(JsonGetInt(part, "Mesh/Coords/NumY") > 0)
    {
        conduit_node_set_path_int(mesh, "coordsets/coords/dims/j", JsonGetInt(part, "Mesh/Coords/NumY"));
        conduit_node_set_path_double(mesh, "coordsets/coords/origin/y", JsonGetDbl(part, "Mesh/Coords/OriginY"));
        conduit_node_set_path_double(mesh, "coordsets/coords/spacing/dy", JsonGetDbl(part, "Mesh/Coords/DeltaY"));
    }

    if(JsonGetInt(part, "Mesh/Coords/NumZ") > 1)
    {
        conduit_node_set_path_int(mesh, "coordsets/coords/dims/k", JsonGetInt(part, "Mesh/Coords/NumZ"));
        conduit_node_set_path_double(mesh, "coordsets/coords/origin/z", JsonGetDbl(part, "Mesh/Coords/OriginZ"));
        conduit_node_set_path_double(mesh, "coordsets/coords/spacing/dz", JsonGetDbl(part, "Mesh/Coords/DeltaZ"));
    }

    /* Topology */
    sprintf(key, "topologies/%s/coordset", topoName);
    conduit_node_set_path_char8_str(mesh, key, "coords");
    sprintf(key, "topologies/%s/type", topoName);
    conduit_node_set_path_char8_str(mesh, key, "uniform");

    /* offset into global -- do we need it? conduit_node_set_path_int(mesh, "topologies/topo/origin/i0", 0); */
}

/*-----------------------------------------------------------------------------*/
static void
json_object_to_blueprint_rectilinear(json_object *part, conduit_node *mesh,
    const char *meshName, const char *topoName)
{
    char key[MAXKEYLEN];
    json_object *json_x = NULL, *json_y = NULL, *json_z = NULL;

    /* Create the coordinates for a rectilinear mesh, zero copy. */
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "rectilinear");
    json_x = JsonGetObj(part, "Mesh/Coords/XAxisCoords");
    if(json_x != NULL)
    {
        double *data = (double *)json_object_extarr_data(json_x);
        if(data != NULL)
        {
            int n = JsonGetInt(part, "Mesh/LogDims", 0);
            conduit_node_set_path_external_double_ptr(mesh, "coordsets/coords/values/x", data, n);
        }
    }
    json_y = JsonGetObj(part, "Mesh/Coords/YAxisCoords");
    if(json_y != NULL)
    {
        double *data = (double *)json_object_extarr_data(json_y);
        if(data != NULL)
        {
            int n = JsonGetInt(part, "Mesh/LogDims", 1);
            conduit_node_set_path_external_double_ptr(mesh, "coordsets/coords/values/y", data, n);
        }
    }
    json_z = JsonGetObj(part, "Mesh/Coords/ZAxisCoords");
    if(json_z != NULL)
    {
        double *data = (double *)json_object_extarr_data(json_z);
        if(data != NULL)
        {
            int n = JsonGetInt(part, "Mesh/LogDims", 2);
            conduit_node_set_path_external_double_ptr(mesh, "coordsets/coords/values/z", data, n);
        }
    }

    /* Topology */
    sprintf(key, "topologies/%s/coordset", topoName);
    conduit_node_set_path_char8_str(mesh, key, "coords");
    sprintf(key, "topologies/%s/type", topoName);
    conduit_node_set_path_char8_str(mesh, key, "rectilinear");
}

/*-----------------------------------------------------------------------------*/
static void
json_object_to_blueprint_curvilinear(json_object *part, conduit_node *mesh,
    const char *meshName, const char *topoName)
{
    char key[MAXKEYLEN];
    const char *ijkNames[] = {"i", "j", "k"};
    const char *jsonAxisNames[] = {"Mesh/Coords/XCoords", "Mesh/Coords/YCoords", "Mesh/Coords/ZCoords"};
    const char *conduitAxisNames[] = {"x", "y", "z", "r", "z", "", "r", "theta", "phi"};
    int i, axisOffset = 0, dims[3] = {1,1,1}, ndims = 1;

    /* Pick the names we'll use for the coordinate axes. */
    if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "X,Y,Z") == 0)
        axisOffset = 0;
    else if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "R,Z") == 0)
        axisOffset = 3;
    else if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "R,Theta,Phi") == 0)
        axisOffset = 6;

    /* Get the dimensions. */
    ndims = JsonGetInt(part, "Mesh/TopoDim");
    for(i = 0; i < ndims; ++i)
    {
        int v = JsonGetInt(part, "Mesh/LogDims", i);
        if(v > 1)
            dims[i] = v;
    }

    /* Create the coordinates for a curvilinear mesh, zero copy. */
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "explicit");
    for(i = 0; i < ndims; ++i)
    {
        json_object *arr = JsonGetObj(part, jsonAxisNames[i]);
        if(arr != NULL)
        {
            double *data = (double *)json_object_extarr_data(arr);
            if(data != NULL)
            {
                int nnodes = dims[0]*dims[1]*dims[2];
                sprintf(key, "coordsets/coords/values/%s", conduitAxisNames[axisOffset + i]);
                conduit_node_set_path_external_double_ptr(mesh, key, data, nnodes);
            }
        }
    }

    /* Topology */
    sprintf(key, "topologies/%s/coordset", topoName);
    conduit_node_set_path_char8_str(mesh, key, "coords");
    sprintf(key, "topologies/%s/type", topoName);
    conduit_node_set_path_char8_str(mesh, key, "structured");
    for(i = 0; i < ndims; ++i)
    {
        sprintf(key, "topologies/%s/elements/dims/%s", topoName, ijkNames[i]);
        conduit_node_set_path_int(mesh, key, dims[i]-1);
    }
}

/*-----------------------------------------------------------------------------*/
static conduit_node *
json_object_to_blueprint_unstructured(json_object *part, const char *meshName)
{
    /* Create the mesh */
    conduit_node *mesh = conduit_node_create();
    return mesh;
}

/*-----------------------------------------------------------------------------*/
static conduit_node *
json_object_to_blueprint_arbitrary(json_object *part, const char *meshName)
{
    /* Create the mesh */
    conduit_node *mesh = conduit_node_create();
    return mesh;
}

/*-----------------------------------------------------------------------------*/
static void
json_object_add_variables_to_conduit_mesh(json_object *part,
    conduit_node *mesh,
    const char *meshName,
    const char *topoName)
{
    char key[MAXKEYLEN];
    json_object *vars_array = JsonGetObj(part, "Vars");
    if(vars_array != NULL)
    {
        /* Compute sizes */
        int i, ndims, ncells = 1, nnodes = 1, dims[3] = {1,1,1}, cdims[3] = {1,1,1};
        ndims = JsonGetInt(part, "Mesh/GeomDim");
        for(i = 0; i < ndims; ++i)
        {
            dims[i] = JsonGetInt(part, "Mesh/LogDims", i);
            cdims[i] = dims[i]-1;
            nnodes *= dims[i];
            ncells *= cdims[i];
        }

        /* Iterate over the fields in the part and add them to the mesh as Conduit fields.*/
        int cent, dtype;
        for (i = 0; i < json_object_array_length(vars_array); i++)
        {
            json_object *varobj = NULL, *dataobj = NULL;
            if((varobj = json_object_array_get_idx(vars_array, i)) != NULL)
            {
                int n = 0;

                /* association */
                sprintf(key, "fields/%s/association", JsonGetStr(varobj, "name"));
                if(strcmp(JsonGetStr(varobj, "centering"),"zone") == 0)
                {
                    conduit_node_set_path_char8_str(mesh, key, "element");
                    n = ncells;
                }
                else
                {
                    conduit_node_set_path_char8_str(mesh, key, "vertex");
                    n = nnodes;
                }
                /* type */
                sprintf(key, "fields/%s/type", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "scalar");
#if 0
                /* volume_dependent */
                sprintf(key, "fields/%s/volume_dependent", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "false"); /* Does MACSio know? */
#endif
                /* topology */
                sprintf(key, "fields/%s/topology", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, topoName);
#if 0
                /* grid_function */
                sprintf(key, "fields/%s/grid_function", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "braid");
#endif
/*CHECK: grid_function */
                /* values */
                sprintf(key, "fields/%s/values", JsonGetStr(varobj, "name"));
                dataobj = JsonGetObj(varobj, "data");
                if(dataobj != NULL)
                {
                    /* Pass data zero copy */
                    json_extarr_type dtype = json_object_extarr_type(dataobj);
                    if(dtype == json_extarr_type_flt32)
                    {
                        conduit_node_set_path_external_float_ptr(mesh, key, 
                            (float *)json_object_extarr_data(dataobj), n);
                    }
                    else if(dtype == json_extarr_type_flt64)
                    {
                        conduit_node_set_path_external_double_ptr(mesh, key, 
                            (double *)json_object_extarr_data(dataobj), n);
                    }
                    else if(dtype == json_extarr_type_int32)
                    {
                        conduit_node_set_path_external_int_ptr(mesh, key, 
                            (int *)json_object_extarr_data(dataobj), n);
                    }
                    else if(dtype == json_extarr_type_int64)
                    {
                        conduit_node_set_path_external_long_ptr(mesh, key, 
                            (long *)json_object_extarr_data(dataobj), n);
                    }
                    else if(dtype == json_extarr_type_byt08)
                    {
                        conduit_node_set_path_external_char_ptr(mesh, key, 
                            (char *)json_object_extarr_data(dataobj), n);
                    }
                    else
                    {
                        MACSIO_LOG_MSG(Warn, ("Unsupported field data type"));
                    }
                }
                else
                {
                    MACSIO_LOG_MSG(Warn, ("Could not get dataobj."));
                }
            }
            else
            {
                MACSIO_LOG_MSG(Warn, ("Could not get varobj."));
            }
        } /* end for */
    }
    else
    {
        MACSIO_LOG_MSG(Warn, ("Could not get Vars."));
    }
}

/*-----------------------------------------------------------------------------*/
/*Q: If there were multiple parts or meshes, would I need to prepend a mesh name to the path? */
static conduit_node *
part_to_conduit(json_object *part, conduit_node *mesh, const char *meshName)
{
    const char *topoName = "mesh";

    if (!strcmp(meshName, "uniform"))
        json_object_to_blueprint_uniform(part, mesh, meshName, topoName);
    else if (!strcmp(meshName, "rectilinear"))
        json_object_to_blueprint_rectilinear(part, mesh, meshName, topoName);
    else if (!strcmp(meshName, "curvilinear"))
        json_object_to_blueprint_curvilinear(part, mesh, meshName, topoName);
#if 0
    else if (!strcmp(meshName, "ucdzoo"))
        mesh = json_object_to_blueprint_unstructured(part, meshName);
    else if (!strcmp(meshName, "arbitrary"))
        mesh = json_object_to_blueprint_arbitrary(part, meshName);
#endif

    json_object_add_variables_to_conduit_mesh(part, mesh, meshName, topoName);

    return mesh;
}

/*-----------------------------------------------------------------------------*/
/* This is an entry point from MACSio to the conduit plugin. */
static void
main_dump(int argi, int argc, char **argv, json_object *main_obj,
    int dumpn, double dumpt)
{
    int i, rank, size;
    json_object *parts = NULL;
    char filename[1024], rootname[1024];
    const char *meshName = "Mesh";

    /* What are we getting for these values? */
    printf("argi=%d\n", argi);
    printf("argc=%d\n", argc);
    printf("argv={");
    for(i = 0; i < argc; ++i)
        printf("%s, ", argv[i]);
    printf("}\n");
    printf("dumpn=%d\n", dumpn);
    printf("dumpt=%lg\n", dumpt);

    /* Now, WTF is main_obj? Args and data? */

    /* MPI rank and size. */
    rank = JsonGetInt(main_obj, "parallel/mpi_rank");
    size = JsonGetInt(main_obj, "parallel/mpi_size");
    printf("rank=%d\nsize=%d\n", rank, size);

    parts = JsonGetObj(main_obj, "problem/parts");
    if(parts != NULL)
    {
        int numParts = 0, nnodes = 0;
        conduit_node **nodes = NULL;

        /* Build Conduit/Blueprint representations of the MACSio json data. */
        numParts = json_object_array_length(parts);
        printf("numParts=%d\n", numParts);

        nodes = (conduit_node **)malloc(numParts * sizeof(conduit_node*));
        if(nodes != NULL)
        {
            const char *meshName = NULL;

            for (i = 0; i < numParts; i++)
            {
                json_object *this_part = JsonGetObj(main_obj, "problem/parts", i);
                if(this_part != NULL)
                {
                    conduit_node *mesh = NULL;
                    if(meshName == NULL)
                        meshName = JsonGetStr(this_part, "Mesh/MeshType");

                    /* Create a node for this mesh part. */
                    nodes[nnodes] = conduit_node_create();

                    /* If we're using a format that will generate blueprint data
                     * then we need the blueprint index file for VisIt to be able
                     * to read the file. For Silo, VisIt would use the Silo reader.
                     * Creating the index via conduit_blueprint_mesh_generate_index
                     * makes a mesh-named set in the index so we'd need to make
                     * sure that the data has that same level. Otherwise the index
                     * and the data don't match and VisIt chokes.
                     */
                    if(strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
                        mesh = conduit_node_fetch(nodes[nnodes], meshName);
                    else
                        mesh = nodes[nnodes];

                    conduit_node_set_path_double(mesh, "state/time", dumpt);
                    conduit_node_set_path_int(mesh, "state/cycle", dumpn);
                    part_to_conduit(this_part, mesh, meshName);

                    nnodes++;
                }
            }

#define SIF 0
            if(SIF)
            {
                /* SIF: Write all of the parts to a single file. */

                /* NOTE: We want this for ADIOS. We also want to represent that 
                         we're adding a new time step to the dataset in hopes 
                         we can have in transit work.
                 */

                /*Q: Do we build up a conduit node to contain the other conduit nodes
                     before passing the data into relay?
                 */
            }
            else
            {
                conduit_node *root = NULL, *bp_idx = NULL, *cindex = NULL, *props = NULL;

                /* MIF: Write all of the parts to separate files. */
                /* The comment about in transit probably applies here too. */
                for(i = 0; i < nnodes; ++i)
                {
                    int ver = 0;
                    conduit_node *info = conduit_node_create();
/*#define DEBUG_PRINT*/
#ifdef DEBUG_PRINT
                    /* Print the data and verify it with Blueprint. */
                    printf("Part %d/%d\n", i+1, nnodes);
                    conduit_node_print(nodes[i]);
                    printf("==================================================================\n");
                    ver = conduit_blueprint_verify("mesh",
                                                   nodes[i],
                                                   info);
                    printf("verify = %d\n", ver);
                    conduit_node_print(info);
                    printf("==================================================================\n");
                    conduit_node_destroy(info);

                    printf("SAVING DATA\n==================================================================\n");
#endif
                    /* Save the data using Conduit. */
                    sprintf(filename, "part.%03d.%03d.%03d%s", rank, i, dumpn, g_preferred_protocol_ext);
                    conduit_relay_io_save2(nodes[i], filename, g_preferred_protocol);
                }

                /* Make the Blueprint index file. */
                if(strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
                {
                    root  = conduit_node_create();
                    bp_idx = conduit_node_fetch(root, "blueprint_index");
                    cindex = conduit_node_fetch(bp_idx, meshName);

#ifdef DEBUG_PRINT
                    printf("GENERATE INDEX\n==================================================================\n");
#endif
                    conduit_blueprint_mesh_generate_index(
                        ((strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
                         ? conduit_node_fetch(nodes[0], meshName)
                         : nodes[0]
                        ),
                        meshName,
                        numParts,
                        cindex);

                    /* We need to store the protocol used to save the data files. */
                    conduit_node_set_path_char8_str(root, "protocol/name", g_preferred_protocol);
                    conduit_node_set_path_char8_str(root, "protocol/version",
                    conduit_node_as_char8_str(conduit_node_fetch(g_conduit_about, "version")));

                    conduit_node_set_path_int(root, "number_of_files", numParts);
                    conduit_node_set_path_int(root, "number_of_trees", 1);
                    conduit_node_set_path_char8_str(root, "file_pattern", filename);
                    conduit_node_set_path_char8_str(root, "tree_pattern", "/");
#ifdef DEBUG_PRINT
                    printf("SAVE INDEX\n==================================================================\n");
#endif
                    /* Save the index using json. */
                    sprintf(rootname, "part.%03d.root", dumpn);
                    conduit_relay_io_save2(root, rootname, "json");
                    conduit_node_destroy(root);
                }
            }

            /* Free nodes */
            for(i = 0; i < nnodes; ++i)
                conduit_node_destroy(nodes[i]);
            free(nodes);
        }
    }

    /* Call relay with the conduit node. */
}

/* Read API. */
static void
main_load(int argi, int argc, char **argv, char const *path, json_object *main_obj,
    json_object **data_read_obj)
{
    int my_rank = JsonGetInt(main_obj, "parallel/mpi_rank");
    int mpi_size = JsonGetInt(main_obj, "parallel/mpi_rank");
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_info(const char *msg, const char *file, int line)
{
    printf("Information:\n  File: %s\n  Line: %d\n  Message: %s\n", file, line, msg);
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_warning(const char *msg, const char *file, int line)
{
    printf("Warning:\n  File: %s\n  Line: %d\n  Message: %s\n", file, line, msg);
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_error(const char *msg, const char *file, int line)
{
    printf("Error:\n  File: %s\n  Line: %d\n  Message: %s\n", file, line, msg);
}

/*-----------------------------------------------------------------------------*/
static int register_this_interface()
{
    MACSIO_IFACE_Handle_t iface;

    if (strlen(iface_name) >= MACSIO_IFACE_MAX_NAME)
        MACSIO_LOG_MSG(Die, ("Interface name \"%s\" too long", iface_name));

//#warning DO Conduit LIB WIDE (DEFAULT) INITIALIZATIONS HERE
    conduit_utils_set_info_handler(macsio_conduit_print_info);
    conduit_utils_set_warning_handler(macsio_conduit_print_warning);
    conduit_utils_set_error_handler(macsio_conduit_print_error);

    /* Get some information about the library. */
    g_conduit_about = conduit_node_create();
    conduit_about(g_conduit_about);
    g_conduit_relay_about = conduit_node_create();
    conduit_relay_about(g_conduit_relay_about);
    set_preferred_protocol();

    /* Populate information about this plugin */
    strcpy(iface.name, iface_name);
    strcpy(iface.ext, iface_ext);
    iface.dumpFunc = main_dump;
    iface.loadFunc = main_load;
    iface.processArgsFunc = process_args;

    /* Register this plugin */
    if (!MACSIO_IFACE_Register(&iface))
        MACSIO_LOG_MSG(Die, ("Failed to register interface \"%s\"", iface_name));

    return 0;
}

/* this one statement is the only statement requiring compilation by
   a C++ compiler. That is because it involves initialization and non
   constant expressions (a function call in this case). This function
   call is guaranteed to occur during *initialization* (that is before
   even 'main' is called) and so will have the effect of populating the
   iface_map array merely by virtue of the fact that this code is linked
   with a main. */
static int dummy = register_this_interface();

/*!@}*/

/*!@}*/
