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

#include <conduit/conduit.h>
#include <conduit/conduit_blueprint.h>
#include <conduit/conduit_utils.h>

#ifdef HAVE_MPI
  #include <mpi.h>
  #include <conduit/conduit_relay_mpi_io.h>
  #define CONDUIT_RELAY_MPI_IO_INITIALIZE(COMM)           conduit_relay_mpi_io_initialize(MPI_Comm_c2f(COMM))
  #define CONDUIT_RELAY_MPI_IO_ABOUT(A,COMM)              conduit_relay_mpi_io_about(A,MPI_Comm_c2f(COMM))
  #define CONDUIT_RELAY_MPI_IO_SAVE2(A,B,C,COMM)          conduit_relay_mpi_io_save2(A,B,C,MPI_Comm_c2f(COMM))
  #define CONDUIT_RELAY_MPI_IO_SAVE3(A,B,C,D,COMM)        conduit_relay_mpi_io_save3(A,B,C,D,MPI_Comm_c2f(COMM))
  #define CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(A,B,C,D,COMM) conduit_relay_mpi_io_save_merged3(A,B,C,D,MPI_Comm_c2f(COMM))
#else
  #include <conduit/conduit_relay.h>
  #define CONDUIT_RELAY_MPI_IO_INITIALIZE(COMM)           conduit_relay_io_initialize()
  #define CONDUIT_RELAY_MPI_IO_ABOUT(A,COMM)              conduit_relay_io_about(A)
  #define CONDUIT_RELAY_MPI_IO_SAVE2(A,B,C,COMM)          conduit_relay_io_save2(A,B,C)
  #define CONDUIT_RELAY_MPI_IO_SAVE3(A,B,C,D,COMM)        conduit_relay_io_save3(A,B,C,D)
  #define CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(A,B,C,D,COMM) conduit_relay_io_save_merged3(A,B,C,D)
#endif

#define MAXKEYLEN 200
#define MAXMSGLEN 2048
#define MAXFILENAMELEN 1024
#define TREE_PATTERN "domain%07d"
#define DIVIDER_STRING "==================================================================\n"

/* Disable debugging messages */
/*#define DEBUG_PRINT*/

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
static char const *iface_ext = "conduit";

/* Static vars */
static conduit_node *g_conduit_about = NULL;
static conduit_node *g_conduit_relay_about = NULL;
static char g_preferred_protocol[100];
static char g_preferred_protocol_ext[100];

/*-----------------------------------------------------------------------------*/
/* We access this one through a function so we can defer its creation until
 * later in the execution process. We store options for conduit protocols
 * into this object. We add the options part of it to meshes that we later
 * write out.
 */
conduit_node *
relay_about(void)
{
    if(g_conduit_relay_about == NULL)
    {
        g_conduit_relay_about = conduit_node_create();
        CONDUIT_RELAY_MPI_IO_ABOUT(g_conduit_relay_about, MACSIO_MAIN_Comm);
    }
    return g_conduit_relay_about;
}

/*-----------------------------------------------------------------------------*/
void
macsio_conduit_finalize(void)
{
    if(g_conduit_about != NULL)
    {
        conduit_node_destroy(g_conduit_about);
        g_conduit_about = NULL;
    }
    if(g_conduit_relay_about != NULL)
    {
        conduit_node_destroy(g_conduit_relay_about);
        g_conduit_relay_about = NULL;
    }
}

/*-----------------------------------------------------------------------------*/
static void
set_preferred_protocol(void)
{
    strcpy(g_preferred_protocol, "conduit_bin");
    strcpy(g_preferred_protocol_ext, "bin");
}

/*-----------------------------------------------------------------------------*/
static char **macsio_conduit_protocols(int *nprotocols)
{
    conduit_node *iop = NULL;
    char **protocols = NULL;
    iop = conduit_node_fetch(relay_about(), "io/protocols");
    if(iop != NULL)
    {
        int i;
        *nprotocols = conduit_node_number_of_children(iop);
        protocols = (char **)malloc(sizeof(char*) * (*nprotocols));
        for(i = 0; i < *nprotocols; ++i)
        {
            char *p = conduit_node_path(conduit_node_child(iop, i));
            if(p != NULL)
            {
                /* p is of the form "io/protocols/conduit_bin" */
                char *start, *slash;
                start = p;
                slash = strrchr(p, '/');
                if(slash != NULL)
                    start = slash + 1;
                protocols[i] = strdup(start);
                free(p);
            }
            else
                protocols[i] = strdup("");
        }
    }
    else
    {  
        *nprotocols = 0;
    }
    return protocols;
}

/*-----------------------------------------------------------------------------*/
static char *macsio_conduit_protocol_help_str(void)
{
    char **protocols = NULL;
    char *protos = NULL;
    int  i, nprotocols = 0, plen = 0;
    const char *cp = "Output protocol for Conduit: [";
    char *c_protocol = g_preferred_protocol;

    protocols = macsio_conduit_protocols(&nprotocols);
    plen += strlen(cp) + 2 + 1;
    for(i = 0; i < nprotocols; ++i)
    {
        plen += strlen(protocols[i]) + 2 + 1;
    }
    protos = (char *)malloc(sizeof(char) * plen);
    if(protos != NULL)
    {
        char *p = protos;
        strcpy(p, cp);
        p += strlen(cp);
        for(i = 0; i < nprotocols; ++i)
        {
            if(i < nprotocols-1)
            {
                sprintf(p, "%s, ", protocols[i]);
                p += strlen(protocols[i]) + 2;
            }
            else
            {
                strcpy(p, protocols[i]);
                p += strlen(protocols[i]);
            }
        }
        strcpy(p, "]\n");
    }

    for(i = 0; i < nprotocols; ++i)
        free(protocols[i]);
    if(nprotocols > 0)
        free(protocols);

    return protos;
}

/*-----------------------------------------------------------------------------*/
static int macsio_conduit_count_options(conduit_node *n)
{
     int i, count = 0, nc = conduit_node_number_of_children(n);
     if(nc > 0)
     {
         for(i = 0; i < nc; ++i)
             count += macsio_conduit_count_options(conduit_node_child(n, i));
     }
     else
     {
         count = 1;
     }

     return count;
}

/*-----------------------------------------------------------------------------*/
static void macsio_conduit_get_options(conduit_node *n, char **opt, int *nopt)
{
     int i, err = 0, count = 0, nc = conduit_node_number_of_children(n);
     char *path, *key;
     if(nc > 0)
     {
         for(i = 0; i < nc; ++i)
             macsio_conduit_get_options(conduit_node_child(n, i), opt, nopt);
     }
     else
     {
         path = conduit_node_path(n);
         key = path + strlen("io/options/");
         opt[*nopt] = strdup(key);        
         free(path);
         *nopt = *nopt + 1;
     }
}

/*-----------------------------------------------------------------------------*/
static void macsio_conduit_add_options(conduit_node *n, char **opt, int nopt, char *p)
{
    int i;
    for(i = 0; i < nopt; ++i)
    {
        int err = 0;
        conduit_node *onode = NULL;
        if((onode = conduit_node_fetch(n, opt[i])) != NULL)
        {
            if(conduit_datatype_is_string(conduit_node_dtype(onode)))
                sprintf(p, "\t%s = \"%s\"\n", opt[i], conduit_node_as_char8_str(onode));
            else if(conduit_datatype_is_int(conduit_node_dtype(onode)))
                sprintf(p, "\t%s = %d\n", opt[i], conduit_node_as_int(onode));
            else if(conduit_datatype_is_float(conduit_node_dtype(onode)))
                sprintf(p, "\t%s = %f\n", opt[i], conduit_node_as_float(onode));
            else if(conduit_datatype_is_double(conduit_node_dtype(onode)))
                sprintf(p, "\t%s = %lg\n", opt[i], conduit_node_as_double(onode));
            else
            {
                /* TODO: fill in other types ...*/
                err = 1;
            }
            if(err == 0)
                p += strlen(p);
        }
    }
}

/*-----------------------------------------------------------------------------*/
static char *macsio_conduit_option_help_str(conduit_node *n, char **opt, int nopt)
{
     const char *co =
         "Pass an option to the Conduit protocol. Options take the form of a path\n"
         "followed by a value. Supported options and defaults are: \n\n";
     char *opts = NULL, *p = NULL;
     int len;

     len = strlen(co)+1 + nopt * 80; /* 80 char per option */
     p = opts = (char *)malloc(sizeof(char) * len);
     strcpy(p, co);
     p += strlen(co);
     macsio_conduit_add_options(n, opt, nopt, p);

     return opts;
}

/*-----------------------------------------------------------------------------*/
static int process_option(conduit_node *n, const char *option, const char *value)
{
    int retval = 0;
    conduit_node *onode = NULL;
    if((onode = conduit_node_fetch(n, option)) != NULL)
    {
        if(conduit_datatype_is_string(conduit_node_dtype(onode)))
            conduit_node_set_path_char8_str(n, option, value);
        else if(conduit_datatype_is_int(conduit_node_dtype(onode)))
        {
            int v;
            if(sscanf(value, "%d", &v) == 1)
                conduit_node_set_path_int(n, option, v);
            else
                retval = 1;
        }
        else if(conduit_datatype_is_float(conduit_node_dtype(onode)))
        {
            float v;
            if(sscanf(value, "%f", &v) == 1)
                conduit_node_set_path_float(n, option, v);
            else
                retval = 1;
        }
        else if(conduit_datatype_is_double(conduit_node_dtype(onode)))
        {
            double v;
            if(sscanf(value, "%lg", &v) == 1)
                conduit_node_set_path_double(n, option, v);
            else
                retval = 1;
        }
        else
        {
            /* TODO: fill in other types ...*/
            retval = 2;
        }
    }
    return retval;
}

/*-----------------------------------------------------------------------------*/
/* This is an entry point from MACSio to the conduit plugin. */
static int process_args(int argi, int argc, char *argv[])
{
    const MACSIO_CLARGS_ArgvFlags_t argFlags = {MACSIO_CLARGS_WARN, MACSIO_CLARGS_TOMEM};
    char *protos_str = NULL, *options_str = NULL;
    char *c_protocol, *c_option;
    char dummy[200];
    char **opt = NULL;
    int i, nopt = 0;

    conduit_node *options = NULL;

    c_protocol = g_preferred_protocol;
    c_option = dummy;
#if 1

#endif
    /* Make sure that we have initialized relay if we have not done it. */
    CONDUIT_RELAY_MPI_IO_INITIALIZE(MACSIO_MAIN_Comm);

    /* Get options. */
    options = conduit_node_fetch(relay_about(), "io/options");
    if(options != NULL)
    {
        nopt = macsio_conduit_count_options(options);
        if(nopt > 0)
        {
            int idx = 0;
            opt = (char **)malloc(sizeof(char*)*(nopt+1));
            macsio_conduit_get_options(options, opt, &idx);
        }
    }

    protos_str = macsio_conduit_protocol_help_str();
    options_str = macsio_conduit_option_help_str(options, opt, nopt);

    /* This thing seems to let MACSio parse args and put values into static global vars in this driver. */
    MACSIO_CLARGS_ProcessCmdline(0, argFlags, argi, argc, argv,
        "--protocol %s", g_preferred_protocol,
            protos_str,
            &c_protocol,
        "--option %s %s", MACSIO_CLARGS_NODEFAULT,
            options_str,
            &c_option, &c_option, /* We don't parse the options using MACSio */
           MACSIO_CLARGS_END_OF_ARGS);

    free(protos_str);
    free(options_str);

#if 0
    printf("BEFORE:\n");
    conduit_node_print(options);
#endif
    /* Handle --option */
    for(i = argi; i < argc; ++i)
    {
        if(strcmp(argv[i], "--option") == 0)
        {
            const char *option, *value;
            int oret = 0;
            option = argv[i+1];
            value = argv[i+2];
            oret = process_option(options, option, value);
            if(oret == 1)
            {
                MACSIO_LOG_MSG(Die, ("Invalid value for option %s", option));
            }
            else if(oret == 2)
            {
                MACSIO_LOG_MSG(Die, ("Unsupported type for option %s", option));
            }
            i += 2;
        }
    }
#if 0
    printf("AFTER:\n");
    conduit_node_print(options);
#endif

    for(i = 0; i < nopt; ++i)
        free(opt[i]);
    free(opt);

    /* Get a file format extension for the protocol -- this really ought to come from Conduit. */
    if(strcmp(g_preferred_protocol, "conduit_bin") == 0)
        strcpy(g_preferred_protocol_ext, "bin");
    else if(strcmp(g_preferred_protocol, "json") == 0)
        strcpy(g_preferred_protocol_ext, "json");
    else if(strcmp(g_preferred_protocol, "hdf5") == 0)
        strcpy(g_preferred_protocol_ext, "h5");
    else if(strcmp(g_preferred_protocol, "conduit_silo") == 0)
        strcpy(g_preferred_protocol_ext, "silo");
    else if(strcmp(g_preferred_protocol, "conduit_silo_mesh") == 0)
        strcpy(g_preferred_protocol_ext, "silo");
    else if(strcmp(g_preferred_protocol, "mpi") == 0)
        strcpy(g_preferred_protocol_ext, "mpi");
    else if(strcmp(g_preferred_protocol, "adios") == 0)
        strcpy(g_preferred_protocol_ext, "bp");
    else
    {
        set_preferred_protocol();
    }

    return 0;
}

/*-----------------------------------------------------------------------------*/
static void
json_object_to_blueprint_origin(json_object *part, conduit_node *mesh,
    const char *topoName)
{
    char key[MAXKEYLEN];
    const char *originName[3] = {"i0", "j0", "k0"};
    int i;
    json_object *origin_obj = JsonGetObj(part, "GlobalLogOrigin");
    for(i = 0; i < json_object_array_length(origin_obj); ++i)
    {
        int val = json_object_get_int(json_object_array_get_idx(origin_obj, i));
        sprintf(key, "topologies/%s/elements/origin/%s", topoName, originName[i]);
        conduit_node_set_path_int(mesh, key, val);
    }
}

/*-----------------------------------------------------------------------------*/
static int
json_object_to_blueprint_uniform(json_object *part, conduit_node *mesh,
    const char *topoName)
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
    json_object_to_blueprint_origin(part, mesh, topoName);

    return 0;
}

/*-----------------------------------------------------------------------------*/
static int
json_object_to_blueprint_rectilinear(json_object *part, conduit_node *mesh,
    const char *topoName)
{
    char key[MAXKEYLEN];
    json_object *json_x = NULL, *json_y = NULL, *json_z = NULL;
    int i, origin[3] = {0,0,0};

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
            origin[0] = JsonGetInt(part, "GlobalLogOrigin", 0);
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
            origin[1] = JsonGetInt(part, "GlobalLogOrigin", 1);
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
            origin[2] = JsonGetInt(part, "GlobalLogOrigin", 2);
        }
    }

    /* Topology */
    sprintf(key, "topologies/%s/coordset", topoName);
    conduit_node_set_path_char8_str(mesh, key, "coords");
    sprintf(key, "topologies/%s/type", topoName);
    conduit_node_set_path_char8_str(mesh, key, "rectilinear");
    json_object_to_blueprint_origin(part, mesh, topoName);

    return 0;
}

/*-----------------------------------------------------------------------------*/
static void
json_object_to_blueprint_explicit_coords(json_object *part, conduit_node *mesh,
    /*out*/int *ndims, /*out*/int dims[3])
{
    char key[MAXKEYLEN];
    const char *jsonAxisNames[] = {"Mesh/Coords/XCoords", "Mesh/Coords/YCoords", "Mesh/Coords/ZCoords"};
    const char *conduitAxisNames[] = {"x", "y", "z", "r", "z", "", "r", "theta", "phi"};
    int i, axisOffset = 0;

    /* Pick the names we'll use for the coordinate axes. */
    if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "X,Y,Z") == 0)
        axisOffset = 0;
    else if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "R,Z") == 0)
        axisOffset = 3;
    else if(strcmp(JsonGetStr(part, "Mesh/Coords/CoordBasis"), "R,Theta,Phi") == 0)
        axisOffset = 6;

    /* Get the dimensions. */
    *ndims = JsonGetInt(part, "Mesh/TopoDim");
    for(i = 0; i < *ndims; ++i)
    {
        int v = JsonGetInt(part, "Mesh/LogDims", i);
        dims[i] = (v > 1) ? v : 1;
    }

    /* Create the coordinates for a curvilinear mesh, zero copy. */
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "explicit");
    for(i = 0; i < *ndims; ++i)
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
}

/*-----------------------------------------------------------------------------*/
static int
json_object_to_blueprint_curvilinear(json_object *part, conduit_node *mesh,
    const char *topoName)
{
    char key[MAXKEYLEN];
    const char *ijkNames[] = {"i", "j", "k"};
    int i, ndims = 1, dims[3] = {1,1,1};

    /* Add the explicit coordinates. */
    json_object_to_blueprint_explicit_coords(part, mesh, &ndims, dims);

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
    json_object_to_blueprint_origin(part, mesh, topoName);

    return 0;
}

/*-----------------------------------------------------------------------------*/
static int
json_object_to_blueprint_unstructured(json_object *part, conduit_node *mesh,
    const char *topoName)
{
    char key[MAXKEYLEN];
    int ndims = 1, dims[3] = {1,1,1}, err = 0;
    json_object *nodelist = NULL;

    /* Topology */
    sprintf(key, "topologies/%s/elements/shape", topoName);
    if(strcmp(JsonGetStr(part, "Mesh/Topology/ElemType"), "Beam2") == 0)
        conduit_node_set_path_char8_str(mesh, key, "line");
    else if(strcmp(JsonGetStr(part, "Mesh/Topology/ElemType"), "Quad4") == 0)
        conduit_node_set_path_char8_str(mesh, key, "quad");
    else if(strcmp(JsonGetStr(part, "Mesh/Topology/ElemType"), "Hex8") == 0)
        conduit_node_set_path_char8_str(mesh, key, "hex");
    else
    {
        MACSIO_LOG_MSG(Warn, ("Unsupported element type."));
        err = 1;
    }
    if(err == 0)
    {
        /* More Topology */
        sprintf(key, "topologies/%s/coordset", topoName);
        conduit_node_set_path_char8_str(mesh, key, "coords");
        sprintf(key, "topologies/%s/type", topoName);
        conduit_node_set_path_char8_str(mesh, key, "unstructured");

        /* Pass the node list zero copy. */
        nodelist = JsonGetObj(part, "Mesh/Topology/Nodelist");
        sprintf(key, "topologies/%s/elements/connectivity", topoName);
        conduit_node_set_path_external_int_ptr(mesh, key, 
            (int *)json_object_extarr_data(nodelist),
            json_object_extarr_nvals(nodelist));

        /* Add the explicit coordinates. */
        json_object_to_blueprint_explicit_coords(part, mesh, &ndims, dims);
    }

    return err;
}

/*-----------------------------------------------------------------------------*/
static int
json_object_to_blueprint_arbitrary(json_object *part, conduit_node *mesh,
    const char *topoName)
{
    char key[MAXKEYLEN];
    int ndims = 1, dims[3] = {1,1,1};
    json_object *nodelist = NULL, *facelist = NULL;

    MACSIO_LOG_MSG(Warn, ("Conduit does not support arbitrary polyhedral "
                          "meshes. Save only faces."));

    /* Add the explicit coordinates. */
    json_object_to_blueprint_explicit_coords(part, mesh, &ndims, dims);

    /* Topology */
    sprintf(key, "topologies/%s/coordset", topoName);
    conduit_node_set_path_char8_str(mesh, key, "coords");
    sprintf(key, "topologies/%s/type", topoName);
    conduit_node_set_path_char8_str(mesh, key, "unstructured");
    sprintf(key, "topologies/%s/elements/shape", topoName);
    if (ndims == 1)
    {
        conduit_node_set_path_char8_str(mesh, key, "line");
        /* Pass the node list zero copy. */
        nodelist = JsonGetObj(part, "Mesh/Topology/Nodelist");
        sprintf(key, "topologies/%s/elements/connectivity", topoName);
        conduit_node_set_path_external_int_ptr(mesh, key, 
            (int *)json_object_extarr_data(nodelist),
            json_object_extarr_nvals(nodelist));
    }
    else if(ndims == 2)
    {
        int *conn = NULL, *dest = NULL;
        int *nl = NULL, *fl = NULL, i, n, nf;
        nodelist = JsonGetObj(part, "Mesh/Topology/Nodelist");
        facelist = JsonGetObj(part, "Mesh/Topology/Facelist");
        n = json_object_extarr_nvals(facelist);
        nf = n / 4;

        /* Use the facelist and nodelist to make quads from edges. */
        conn = (int *)malloc(n * sizeof(int));
        nl = (int *)json_object_extarr_data(nodelist);
        fl = (int *)json_object_extarr_data(facelist);
        dest = conn;
        for(i = 0; i < nf; ++i)
        {
            int *e0, *e1;
            /* The first 2 edges that make up the shape. */
            e0 = nl + 2*fl[i*4];
            e1 = nl + 2*fl[i*4+1];

            dest[0] = e0[0];
            dest[1] = e1[0];
            dest[2] = e1[1];
            dest[3] = e0[1];

            dest += 4;
        }

        /* Pass the node list but we have to pass a reordered copy. */       
        conduit_node_set_path_char8_str(mesh, key, "quad");
        sprintf(key, "topologies/%s/elements/connectivity", topoName);
        conduit_node_set_path_int_ptr(mesh, key, conn, n);
    }
    else if(ndims == 3)
    {
        conduit_node_set_path_char8_str(mesh, key, "quad");
        /* Pass the node list zero copy. */
        nodelist = JsonGetObj(part, "Mesh/Topology/Nodelist");
        sprintf(key, "topologies/%s/elements/connectivity", topoName);
        conduit_node_set_path_external_int_ptr(mesh, key, 
            (int *)json_object_extarr_data(nodelist),
            json_object_extarr_nvals(nodelist));
    }

    return 0;
}

/*-----------------------------------------------------------------------------*/
static void
json_object_add_variables_to_blueprint_mesh(json_object *part,
    conduit_node *mesh,
    const char *meshName,
    const char *topoName)
{
    char key[MAXKEYLEN];
    json_object *vars_array = JsonGetObj(part, "Vars");
    if(vars_array != NULL)
    {
        int i;
        /* Iterate over the fields in the part and add them to the mesh as Conduit fields.*/
        for (i = 0; i < json_object_array_length(vars_array); i++)
        {
            json_object *varobj = NULL, *dataobj = NULL;
            if((varobj = json_object_array_get_idx(vars_array, i)) != NULL)
            {
                /* association */
                sprintf(key, "fields/%s/association", JsonGetStr(varobj, "name"));
                if(strcmp(JsonGetStr(varobj, "centering"),"zone") == 0)
                    conduit_node_set_path_char8_str(mesh, key, "element");
                else
                    conduit_node_set_path_char8_str(mesh, key, "vertex");

                /* type */
                sprintf(key, "fields/%s/type", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "scalar");

                /* topology */
                sprintf(key, "fields/%s/topology", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, topoName);
#if 0
                /*NOTE: these fields seem optional. */
                /* volume_dependent */
                sprintf(key, "fields/%s/volume_dependent", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "false"); /* Does MACSio know? */

                /* grid_function */
                sprintf(key, "fields/%s/grid_function", JsonGetStr(varobj, "name"));
                conduit_node_set_path_char8_str(mesh, key, "braid");
#endif

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
                            (float *)json_object_extarr_data(dataobj),
                            json_object_extarr_nvals(dataobj));
                    }
                    else if(dtype == json_extarr_type_flt64)
                    {
                        conduit_node_set_path_external_double_ptr(mesh, key, 
                            (double *)json_object_extarr_data(dataobj),
                            json_object_extarr_nvals(dataobj));
                    }
                    else if(dtype == json_extarr_type_int32)
                    {
                        conduit_node_set_path_external_int_ptr(mesh, key, 
                            (int *)json_object_extarr_data(dataobj),
                            json_object_extarr_nvals(dataobj));
                    }
                    else if(dtype == json_extarr_type_int64)
                    {
                        conduit_node_set_path_external_long_ptr(mesh, key, 
                            (long *)json_object_extarr_data(dataobj),
                            json_object_extarr_nvals(dataobj));
                    }
                    else if(dtype == json_extarr_type_byt08)
                    {
                        conduit_node_set_path_external_char_ptr(mesh, key, 
                            (char *)json_object_extarr_data(dataobj),
                            json_object_extarr_nvals(dataobj));
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
static void
json_object_to_blueprint(json_object *part, conduit_node *mesh, const char *meshName)
{
    const char *topoName = "mesh";
    int err = 0;

    if (!strcmp(meshName, "uniform"))
        err = json_object_to_blueprint_uniform(part, mesh, topoName);
    else if (!strcmp(meshName, "rectilinear"))
        err = json_object_to_blueprint_rectilinear(part, mesh, topoName);
    else if (!strcmp(meshName, "curvilinear"))
        err = json_object_to_blueprint_curvilinear(part, mesh, topoName);
    else if (!strcmp(meshName, "ucdzoo"))
        err = json_object_to_blueprint_unstructured(part, mesh, topoName);
    else if (!strcmp(meshName, "arbitrary"))
        err = json_object_to_blueprint_arbitrary(part, mesh, topoName);

    if(err == 0)
        json_object_add_variables_to_blueprint_mesh(part, mesh, meshName, topoName);
}

/*-----------------------------------------------------------------------------*/
#ifdef DEBUG_PRINT
static void
verify_mesh(conduit_node *mesh)
{
    int ver = 0;
    conduit_node *info;
    info = conduit_node_create();
    conduit_node_print(mesh);

    printf(DIVIDER_STRING);
    ver = conduit_blueprint_verify("mesh", mesh, info);

    printf("verify = %d\n", ver);
    conduit_node_print(info);
    printf(DIVIDER_STRING);
    conduit_node_destroy(info);
}
#endif

/*-----------------------------------------------------------------------------*/
int
conduit_relay_io_supports_collective(const char *protocol)
{
    /* Right now, ADIOS is the only collective writer in Conduit. */
    if(strcmp(protocol, "adios") == 0)
        return 1;
    return 0;
}

/*-----------------------------------------------------------------------------*/
/* @brief Get the I/O options for the preferred protocol.
 */
static conduit_node *macsio_conduit_get_protocol_options(void)
{
    char key[MAXKEYLEN];
    conduit_node *p_options = NULL;

    /* Try and get options for the preferred protocol.
     */
    sprintf(key, "io/options/%s", g_preferred_protocol);
    if(conduit_node_has_path(relay_about(), key))
        p_options = conduit_node_fetch(relay_about(), key);
    return p_options;
}

/*-----------------------------------------------------------------------------*/
/* Given a conduit node that contains multiple branches named "domain0", "domain1",
 * ... save the domains to a file through Conduit.
 *
 * MACSio is allowed to pass multiple parts inside the mesh_root node but for
 * many applications, there will be one part per processor.
 */
static void
macsio_conduit_sif_write(json_object *main_obj,
    int dumpn,
    const char *meshName, 
    conduit_node *mesh_root, conduit_node **nodes, int nnodes, int totalNumParts)
{
    int msg, tag = 12345;
    const char *tree = "/";
    char filename[MAXFILENAMELEN], rootname[MAXFILENAMELEN];
    int rank, size;
    const char *ext = NULL;
#ifdef HAVE_MPI
    MPI_Comm splitcomm;
#else
    int splitcomm = 0;
#endif

    rank = JsonGetInt(main_obj, "parallel/mpi_rank");
    size = JsonGetInt(main_obj, "parallel/mpi_size");
    ext = json_object_path_get_string(main_obj, "clargs/fileext");
    if(strcmp(ext, iface_ext) == 0)
        ext = g_preferred_protocol_ext;

#ifdef HAVE_MPI
    /* This comm is split so each split comm has 1 rank. 
     * We use this to write the index file.
     */
    MPI_Comm_split(MACSIO_MAIN_Comm, rank, 0, &splitcomm);
#endif

#ifdef DEBUG_PRINT
    printf("SIF: mesh\n" DIVIDER_STRING);
    conduit_node_print(mesh_root);
#endif

    /* Save the data using Conduit. */
    sprintf(filename, "%s_conduit_%03d.%s",
            json_object_path_get_string(main_obj, "clargs/filebase"),
            dumpn,
            ext);

    if(conduit_relay_io_supports_collective(g_preferred_protocol))
    {
#ifdef DEBUG_PRINT
        printf("SIF: size=%d, nnodes=%d, numParts=%d\n", size, nnodes, totalNumParts);
#endif
        /* The format supports collective writes.
         * Each rank will write to the same file at the same time.
         * NOTE: We write using nodes[] because this means that domains
         *       will not have any domainXXXXX prefix.
         */
        CONDUIT_RELAY_MPI_IO_SAVE3(nodes[0], filename, g_preferred_protocol,
            macsio_conduit_get_protocol_options(), MACSIO_MAIN_Comm);
        /* We add multiple parts at the end, appending to the file. */
        for(int i = 1; i < nnodes; ++i)
        {
            CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(nodes[i], filename, g_preferred_protocol,
                macsio_conduit_get_protocol_options(), MACSIO_MAIN_Comm);
        }
    }
    else
    {
        /* The format supports only serial writes.
         * The ranks will sequentially append to the file. BAD!
         * NOTE: The data are written using prefixes like domainXXXXX so we
         *       can accommodate multiple domains in a single file for
         *       non-collective formats.
         */
        tree = TREE_PATTERN;

#ifdef HAVE_MPI
        /* NOTE: this is our own baton passing. 
         * All ranks contribute to the same file.
         */
        if(rank == 0)
        {
#ifdef DEBUG_PRINT
            printf("CONDUIT_RELAY_MPI_IO_SAVE3 %s on rank %d\n", filename, rank);
#endif
            CONDUIT_RELAY_MPI_IO_SAVE3(mesh_root, filename, g_preferred_protocol,
                macsio_conduit_get_protocol_options(), splitcomm);
            msg = 0;
            MPI_Send(&msg, 1, MPI_INT, rank+1, tag, MACSIO_MAIN_Comm);
        }
        else
        {
            MPI_Status status;
            MPI_Recv(&msg, 1, MPI_INT, rank-1, tag, MACSIO_MAIN_Comm, &status);

#ifdef DEBUG_PRINT
            printf("CONDUIT_RELAY_MPI_IO_SAVE_MERGED3 %s on rank %d\n", filename, rank);
#endif
            /* Let this rank append to the file. */
            CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(mesh_root, filename, g_preferred_protocol,
                macsio_conduit_get_protocol_options(), splitcomm);

            if(rank < size-1)
                MPI_Send(&msg, 1, MPI_INT, rank+1, tag, MACSIO_MAIN_Comm);
        }
#else
        /* no HAVE_MPI so there is only one rank. */
        if(size > 1)
        {
            MACSIO_LOG_MSG(Die, ("MPI is required for Conduit's SIF mode"));
        }
        CONDUIT_RELAY_MPI_IO_SAVE3(mesh_root, filename, g_preferred_protocol,
            macsio_conduit_get_protocol_options(), MACSIO_MAIN_Comm);
#endif
    }
    /*MACSIO_UTILS_RecordOutputFiles(dumpn, filename);*/

    /* Make the index file for protocols that don't make Silo mesh files. */
    if(rank == 0 && strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
    {
        conduit_node *bp_root = NULL, *bp_idx = NULL, *cindex = NULL, *props = NULL;
        bp_root  = conduit_node_create();
        bp_idx = conduit_node_fetch(bp_root, "blueprint_index");
        cindex = conduit_node_fetch(bp_idx, meshName);
#ifdef DEBUG_PRINT
        printf("GENERATE INDEX\n" DIVIDER_STRING);
#endif
        conduit_blueprint_mesh_generate_index(
            ((strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
             ? conduit_node_fetch(nodes[0], meshName)
             : nodes[0]
            ),
            meshName,
            totalNumParts,
            cindex);
        /* We need to store the protocol used to save the data files. */
        conduit_node_set_path_char8_str(bp_root, "protocol/name", g_preferred_protocol);
        conduit_node_set_path_char8_str(bp_root, "protocol/version",
        conduit_node_as_char8_str(conduit_node_fetch(g_conduit_about, "version")));

        conduit_node_set_path_int(bp_root, "number_of_files", 1);
        conduit_node_set_path_int(bp_root, "number_of_trees", totalNumParts);

        conduit_node_set_path_char8_str(bp_root, "file_pattern", filename);
        conduit_node_set_path_char8_str(bp_root, "tree_pattern", tree);
#ifdef DEBUG_PRINT
        printf("SAVE INDEX\n" DIVIDER_STRING);
#endif
        /* Save the index using json. */
        sprintf(rootname, "%s_conduit_%03d.root",
            json_object_path_get_string(main_obj, "clargs/filebase"), 
            dumpn);
        CONDUIT_RELAY_MPI_IO_SAVE2(bp_root, rootname, "json", splitcomm);
        conduit_node_destroy(bp_root);
        /*MACSIO_UTILS_RecordOutputFiles(dumpn, rootname);*/
    }

#ifdef HAVE_MPI
    MPI_Comm_free(&splitcomm);
#endif
}

typedef struct
{
#ifdef HAVE_MPI
    MPI_Comm comm;
#else
    int comm;
#endif
    conduit_node *mesh_root;
} file_data;

static void *CreateConduitFile(const char *fname, const char *nsname, void *userData)
{
    file_data *filedata = (file_data *)userData;

    /* We're the first to have the baton for the file. save. */
    CONDUIT_RELAY_MPI_IO_SAVE3(filedata->mesh_root, fname, g_preferred_protocol,
       macsio_conduit_get_protocol_options(), filedata->comm);

    int *retval = (int *)malloc(sizeof(int));
    *retval = 1; /* Create */
    return (void *) retval;
}

static void *OpenConduitFile(const char *fname, const char *nsname,
                             MACSIO_MIF_ioFlags_t ioFlags, void *userData)
{
    file_data *filedata = (file_data *)userData;

    /* We're next to have the baton for the file. save merged. */
    CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(filedata->mesh_root, fname, g_preferred_protocol,
       macsio_conduit_get_protocol_options(), filedata->comm);

    int *retval = (int *)malloc(sizeof(int));
    *retval = 2; /* Open */
    return (void *) retval;
}

static void CloseConduitFile(void *file, void *userData)
{
    free(file);
}

/*-----------------------------------------------------------------------------*/
/* Given a list of conduit nodes that represent mesh parts (domains), write 
 * each one to a separate file sequentially through Conduit.
 */
static void
macsio_conduit_mif_write(json_object *main_obj,
    int dumpn,
    const char *meshName, 
    conduit_node *mesh_root, conduit_node **nodes, int nnodes, int numParts, int nFiles)
{
    int i, domain_id, np, nt, write_one_file_per_part = 0;
    char filename[MAXFILENAMELEN], rootname[MAXFILENAMELEN];
    const char *tree = NULL;
#ifdef HAVE_MPI
    MPI_Comm splitcomm;
#else
    int splitcomm = 0;
#endif
    int rank, size;
    const char *ext = NULL;

    rank = JsonGetInt(main_obj, "parallel/mpi_rank");
    size = JsonGetInt(main_obj, "parallel/mpi_size");
    ext = json_object_path_get_string(main_obj, "clargs/fileext");
    if(strcmp(ext, iface_ext) == 0)
        ext = g_preferred_protocol_ext;

#ifdef HAVE_MPI
    /* This comm is split so each split comm has 1 rank. */
    MPI_Comm_split(MACSIO_MAIN_Comm, MACSIO_MAIN_Rank, 0, &splitcomm);
#endif

    if(size == 1)
    {
        if(nFiles == numParts)
            write_one_file_per_part = 1;
    }
    else
    {
        write_one_file_per_part = numParts == nFiles;
    }
#if 1
    /* Force 1 file per part for Silo Mesh. */
    if(strcmp(g_preferred_protocol, "conduit_silo_mesh") == 0)
        write_one_file_per_part = 1;
#endif
    if(write_one_file_per_part)
    {
#ifdef DEBUG_PRINT
        printf("MIF: Write one file per part (nnodes=%d, numParts=%d, nFiles=%d)\n", nnodes, numParts, nFiles);
#endif
        /* Each process will write its own files using splitcomm. */
        for(i = 0; i < nnodes; ++i)
        {
            char key[MAXKEYLEN];
            if(strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
                sprintf(key, "%s/state/domain_id", meshName);
            else
                strcpy(key, "state/domain_id");
            domain_id = conduit_node_as_int(conduit_node_fetch(nodes[i], key));
            sprintf(filename, "%s_conduit_%05d_%03d.%s",
                    json_object_path_get_string(main_obj, "clargs/filebase"),
                    domain_id,
                    dumpn,
                    ext);

#ifdef DEBUG_PRINT
            printf("MIF: mesh %d\n" DIVIDER_STRING, domain_id);
            conduit_node_print(nodes[i]);
#endif

            CONDUIT_RELAY_MPI_IO_SAVE3(nodes[i], filename, g_preferred_protocol,
                macsio_conduit_get_protocol_options(), splitcomm);

            /*printf("Saved %s\n", filename);*/
            /*MACSIO_UTILS_RecordOutputFiles(dumpn, filename);*/
        }
        /* We'll have numParts (global number of parts) files*/
        np = numParts;
        nt = 1;
        tree = "/";
    }
    else
    {
#ifdef DEBUG_PRINT
        printf("MIF: Write %d parts to %d files (nnodes=%d)\n", numParts, nFiles, nnodes);
#endif
        /* NOTE: the output looks right for this. VisIt does not read it when 
         * more than one file is produced. Looks like an issue in VisIt.
         */
        if(nFiles > size)
        {
            MACSIO_LOG_MSG(Warn, ("The output number of files %d is greater "
               "than %d ranks. The number of files will be %d", nFiles, size, size));
        }

#ifdef DEBUG_PRINT
        printf("MIF: mesh\n" DIVIDER_STRING);
        conduit_node_print(nodes[i]);
#endif

        /* Create a baton. We do it up here since it looks like low overhead 
         * and lets us compute the group number of the current MPI rank, 
         * which we use for making file names, etc.
         */
        file_data filedata;
        filedata.comm = splitcomm;
        filedata.mesh_root = mesh_root;
        int *file_ptr = NULL;
        MACSIO_MIF_ioFlags_t ioFlags = {MACSIO_MIF_WRITE,
                JsonGetInt(main_obj, "clargs/exercise_scr")&0x1};
        MACSIO_MIF_baton_t *bat = MACSIO_MIF_Init(nFiles, ioFlags, MACSIO_MAIN_Comm, 3,
             CreateConduitFile, OpenConduitFile, CloseConduitFile, &filedata);

        /* Come up with the filename for the group. */
        int group_rank = MACSIO_MIF_RankOfGroup(bat, rank);
        sprintf(filename, "%s_conduit_%05d_%03d.%s",
                json_object_path_get_string(main_obj, "clargs/filebase"),
                group_rank,
                dumpn,
                ext);

        if(conduit_relay_io_supports_collective(g_preferred_protocol))
        {
#ifdef DEBUG_PRINT
            printf("MIF: Collective group write 1 part to %d files.: rank=%d, group_rank=%d\n", nFiles, rank, group_rank);
#endif
            /* We have a collective format so we have groups of MPI ranks
             * that will all want to write to the file at once.
             */
#ifdef HAVE_MPI
            MPI_Comm groupcomm;
            MPI_Comm_split(MACSIO_MAIN_Comm, group_rank, 0, &groupcomm);
#else
            int groupcom = 0;
#endif
            CONDUIT_RELAY_MPI_IO_SAVE3(nodes[0], filename, g_preferred_protocol,
                macsio_conduit_get_protocol_options(), groupcomm);
            /* Append other parts (assumes nnodes is same on all ranks). */
            for(int i = 1; i < nnodes; ++i)
            {
                CONDUIT_RELAY_MPI_IO_SAVE_MERGED3(nodes[i], filename, g_preferred_protocol,
                    macsio_conduit_get_protocol_options(), groupcomm);
            }

#ifdef HAVE_MPI
            MPI_Comm_free(&groupcomm);
#endif
            /* We'll have nFiles files since multiple ranks might write same file. */
            np = nFiles;
            nt = numParts;
            tree = TREE_PATTERN;
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("MIF: Baton-based appends\n");
#endif
            /* The format is not collective so we must append to it serially.
             * Use the split comm for this.
             */

            /* Wait for our turn writing to the file. */
            file_ptr = (int *) MACSIO_MIF_WaitForBaton(bat, filename, 0);

            /* Hand off the baton to the next processor. This winds up closing
             * the file so that the next processor that opens it can be assured
             * of getting a consistent and up to date view of the file's contents. */
            MACSIO_MIF_HandOffBaton(bat, file_ptr);

            /*MACSIO_UTILS_RecordOutputFiles(dumpn, filename);*/

            /* We'll have nFiles files since multiple ranks might write same file. */
            np = (nFiles > size) ? size : nFiles;
            nt = numParts;
            tree = TREE_PATTERN;
        }

        /* We're done using MACSIO_MIF, so finish it off */
        MACSIO_MIF_Finish(bat);
    }

    /* Make the Blueprint index file. */
    if(rank == 0 && strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
    {
        conduit_node *root = NULL, *bp_idx = NULL, *cindex = NULL, *props = NULL;
        root  = conduit_node_create();
        bp_idx = conduit_node_fetch(root, "blueprint_index");
        cindex = conduit_node_fetch(bp_idx, meshName);

#ifdef DEBUG_PRINT
        printf("GENERATE INDEX\n" DIVIDER_STRING);
#endif
        conduit_blueprint_mesh_generate_index(
            ((strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
             ? conduit_node_fetch(nodes[0], meshName)
             : nodes[0]
            ),
            meshName,
            numParts,
            cindex);

        /* Store the protocol used to save the data files. */
        conduit_node_set_path_char8_str(root, "protocol/name", g_preferred_protocol);
        conduit_node_set_path_char8_str(root, "protocol/version",
        conduit_node_as_char8_str(conduit_node_fetch(g_conduit_about, "version")));

        conduit_node_set_path_int(root, "number_of_files", np);
        conduit_node_set_path_int(root, "number_of_trees", nt);

        sprintf(filename, "%s_conduit_%s_%03d.%s",
            json_object_path_get_string(main_obj, "clargs/filebase"),
            ((nFiles > 1) || (write_one_file_per_part && (numParts > 1))) ? "%05d" : "00000",
            dumpn,
            ext);
        conduit_node_set_path_char8_str(root, "file_pattern", filename);
        conduit_node_set_path_char8_str(root, "tree_pattern", tree);

        /* Save the index using json. */
        sprintf(rootname, "%s_conduit_%03d.root",
                json_object_path_get_string(main_obj, "clargs/filebase"),
                dumpn);
#ifdef DEBUG_PRINT
        printf("SAVE INDEX %s\n" DIVIDER_STRING, rootname);
#endif
        CONDUIT_RELAY_MPI_IO_SAVE2(root, rootname, "json", splitcomm);
        conduit_node_destroy(root);
        /*MACSIO_UTILS_RecordOutputFiles(dumpn, rootname);*/
    }

#ifdef HAVE_MPI
    MPI_Comm_free(&splitcomm);
#endif
}

/*-----------------------------------------------------------------------------*/
typedef enum
{
    FileWrite_SIF,
    FileWrite_MIF
} FileWrite_t;

static void
conduit_macsio_get_file_modes(json_object *main_obj, FileWrite_t *writeT, int *nFiles)
{
    json_object *parfmode_obj = json_object_path_get_array(main_obj, "clargs/parallel_file_mode");
    if (parfmode_obj)
    {
        json_object *modestr = json_object_array_get_idx(parfmode_obj, 0);
        json_object *filecnt = json_object_array_get_idx(parfmode_obj, 1);
        if(strcmp(json_object_get_string(modestr), "SIF") == 0)
        {
            *writeT = FileWrite_SIF;
            *nFiles = 1;
        }
        else if(strcmp(json_object_get_string(modestr), "MIF") == 0)
        {
            *writeT = FileWrite_MIF;
            *nFiles = json_object_get_int(filecnt);
        }
        else if(strcmp(json_object_get_string(modestr), "MIFMAX") == 0)
        {
            *writeT = FileWrite_MIF;
            *nFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
        }
    }
    else
    {
        char const * modestr = json_object_path_get_string(main_obj, "clargs/parallel_file_mode");
        if (strcmp(modestr, "SIF") == 0)
        {
            float avg_num_parts = json_object_path_get_double(main_obj, "clargs/avg_num_parts");
            if (avg_num_parts == (float ((int) avg_num_parts)))
            {
                *writeT = FileWrite_SIF;
                *nFiles = 1;
            }
            else
            {
                MACSIO_LOG_MSG(Die, ("Conduit plugin cannot currently handle SIF mode where "
                    "there are different numbers of parts on each MPI rank. "
                    "Set --avg_num_parts to an integral value." ));
            }
        }
        else if (strcmp(modestr, "MIF") == 0 ||
                 strcmp(modestr, "MIFMAX") == 0)
        {
            *writeT = FileWrite_MIF;
            *nFiles = json_object_path_get_int(main_obj, "parallel/mpi_size");
        }
    }
#if 0
    if(*writeT == FileWrite_MIF)
        printf("MIF, %d\n", *nFiles);
    else if(*writeT == FileWrite_MIF)
        printf("SIF, %d\n", *nFiles);
#endif
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
/*printf("main_dump: 0\n");*/
    /* process cl args */
    process_args(argi, argc, argv);
/*printf("main_dump: 1\n");*/

    /* MPI rank and size. */
    rank = JsonGetInt(main_obj, "parallel/mpi_rank");
    size = JsonGetInt(main_obj, "parallel/mpi_size");
/*printf("main_dump: 2\n");*/

    parts = JsonGetObj(main_obj, "problem/parts");
    if(parts != NULL)
    {
        int numParts = 0, totalNumParts = 0, nnodes = 0;
        conduit_node **nodes = NULL;
/*printf("main_dump: 3\n");*/

        /* Build Conduit/Blueprint representations of the MACSio json data. */
        numParts = json_object_array_length(parts);
        nodes = (conduit_node **)malloc(numParts * sizeof(conduit_node*));
        if(nodes != NULL)
        {
            const char *meshName = NULL, *ext = NULL;
            FileWrite_t writeT = FileWrite_MIF;
            int nFiles = size;
/*printf("main_dump: 4\n");*/

#define COMBINE_ROOT
#ifdef COMBINE_ROOT
            /* Make a single mesh_root with domainXX branches under it to contain
             * the domains.
             */
            char domainpath[100];
            conduit_node *mesh_root = NULL;
            mesh_root = conduit_node_create();
            conduit_node_set_path_double(mesh_root, "time", dumpt);
            conduit_node_set_path_int(mesh_root, "cycle", dumpn);
#endif
            for (i = 0; i < numParts; i++)
            {
                json_object *this_part = JsonGetObj(main_obj, "problem/parts", i);
                if(this_part != NULL)
                {
                    conduit_node *mesh = NULL;
                    if(meshName == NULL)
                        meshName = JsonGetStr(this_part, "Mesh/MeshType");

                    /* Create a node for this mesh part. */
#ifdef COMBINE_ROOT
                    sprintf(domainpath, TREE_PATTERN, JsonGetInt(this_part, "Mesh/ChunkID"));
                    nodes[nnodes] = conduit_node_fetch(mesh_root, domainpath);
#else
                    nodes[nnodes] = conduit_node_create();
#endif

                    /* If we're using a format that will generate blueprint data
                     * then we need the blueprint index file for VisIt to be able
                     * to read the file. For "conduit_silo_mesh", VisIt would use 
                     * the Silo reader. Creating the index via conduit_blueprint_mesh_generate_index
                     * makes a mesh-named set in the index so we'd need to make
                     * sure that the data has that same level. Otherwise the index
                     * and the data don't match and VisIt chokes.
                     */
                    if(strcmp(g_preferred_protocol, "conduit_silo_mesh") != 0)
                        mesh = conduit_node_fetch(nodes[nnodes], meshName);
                    else
                        mesh = nodes[nnodes];

                    /* Add cycle, time, domain */
                    conduit_node_set_path_double(mesh, "state/time", dumpt);
                    conduit_node_set_path_int(mesh, "state/cycle", dumpn);
                    conduit_node_set_path_int(mesh, "state/domain_id", 
                        JsonGetInt(this_part, "Mesh/ChunkID"));

                    /* Add the blueprint mesh data. */
                    json_object_to_blueprint(this_part, mesh, meshName);
#if 0
                    if(rank == 0)
                    {
                        conduit_node_print(mesh_root);
                        conduit_node_print(nodes[nnodes]);
                    }
#endif
#ifdef DEBUG_PRINT
                    /* Verify the mesh data. */
                    verify_mesh(mesh);
#endif
                    nnodes++;
                }
            }
/*printf("main_dump: 5\n");*/

            /* Get the file extension. If it has been overridden by the user,
             * then use that. Otherwise, use the file extension for our
             * preferred protocol.
             */
            ext = json_object_path_get_string(main_obj, "clargs/fileext");
            if(strcmp(ext, iface_ext) == 0)
                ext = g_preferred_protocol_ext;
/*printf("main_dump: 6\n");*/

#ifdef HAVE_MPI
            /* Sum the total number of parts across ranks. */
            totalNumParts = numParts;
            MPI_Allreduce(&numParts, &totalNumParts, 1, 
                          MPI_INT, MPI_SUM, MACSIO_MAIN_Comm);
#else
            totalNumParts = numParts;
#endif
/*printf("main_dump: 7\n");*/

#ifdef DEBUG_PRINT
            printf("SAVING DATA\n" DIVIDER_STRING);
#endif
            /* Get the file mode and number of files to make. */
            conduit_macsio_get_file_modes(main_obj, &writeT, &nFiles);
            if(writeT == FileWrite_SIF)
            {
                /* SIF: Write all of the parts to a single file. */
                macsio_conduit_sif_write(main_obj, dumpn,
                    meshName,
                    mesh_root, nodes, nnodes, totalNumParts);
            }
            else
            {
                /* MIF: Write all of the parts to some number of separate files. */
                macsio_conduit_mif_write(main_obj, dumpn,
                    meshName,
                    mesh_root, nodes, nnodes, totalNumParts, nFiles);
            }
/*printf("main_dump: 8\n");*/

#ifdef COMBINE_ROOT
            conduit_node_destroy(mesh_root);
#else
            /* Free nodes */
            for(i = 0; i < nnodes; ++i)
                conduit_node_destroy(nodes[i]);
#endif
            free(nodes);
        }
    }
/*printf("main_dump: 9\n");*/

    /* Call relay with the conduit node. */
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_info(const char *msg, const char *file, int line)
{
    char buf[MAXMSGLEN];
    snprintf(buf, MAXMSGLEN, "File: %s\n  Line: %d\n  Information: %s", file, line, msg);
    printf("%s\n", buf);
    MACSIO_LOG_MSG(Info, (buf));
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_warning(const char *msg, const char *file, int line)
{
    char buf[MAXMSGLEN];
    snprintf(buf, MAXMSGLEN, "File: %s\n  Line: %d\n  Warning: %s\n", file, line, msg);
    printf("%s\n", buf);
    MACSIO_LOG_MSG(Warn, (buf));
}

/*-----------------------------------------------------------------------------*/
static void 
macsio_conduit_print_error(const char *msg, const char *file, int line)
{
    char buf[MAXMSGLEN];
    snprintf(buf, MAXMSGLEN, "File: %s\n  Line: %d\n  Error: %s\n", file, line, msg);
    printf("%s\n", buf);
    MACSIO_LOG_MSG(Err, (buf));
}

/*-----------------------------------------------------------------------------*/
static int register_this_interface()
{
    MACSIO_IFACE_Handle_t iface;

    if (strlen(iface_name) >= MACSIO_IFACE_MAX_NAME)
        MACSIO_LOG_MSG(Die, ("Interface name \"%s\" too long", iface_name));

//#warning DO Conduit LIB WIDE (DEFAULT) INITIALIZATIONS HERE
    /* Install some message handlers for Conduit so we can print messages. */
    conduit_utils_set_info_handler(macsio_conduit_print_info);
    conduit_utils_set_warning_handler(macsio_conduit_print_warning);
    conduit_utils_set_error_handler(macsio_conduit_print_error);

    /* Get some information about the library. */
    g_conduit_about = conduit_node_create();
    conduit_about(g_conduit_about);
    set_preferred_protocol();
    atexit(macsio_conduit_finalize);

    /* Populate information about this plugin */
    strcpy(iface.name, iface_name);
    strcpy(iface.ext, iface_ext);
    iface.dumpFunc = main_dump;
    iface.processArgsFunc = process_args;

    /* Register this plugin */
    if (!MACSIO_IFACE_Register(&iface))
        MACSIO_LOG_MSG(Die, ("Failed to register interface \"%s\"", iface_name));

    return 0;
}

#if 0
/* this one statement is the only statement requiring compilation by
   a C++ compiler. That is because it involves initialization and non
   constant expressions (a function call in this case). This function
   call is guaranteed to occur during *initialization* (that is before
   even 'main' is called) and so will have the effect of populating the
   iface_map array merely by virtue of the fact that this code is linked
   with a main. */
static int dummy = register_this_interface();
#else
/* Provide a function we can call from macsio_main once MPI is initialized. */
void register_conduit_interface(void)
{
    register_this_interface();
}
#endif

/*!@}*/

/*!@}*/
