#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include "common.h"
#include "hashtable.h"

int main (void){
    char path [] = {"a06161.1.flds.zfem"};
    int npoin,nelem,*elems;
    double *ptxyz;
    read_zfem(path,&npoin,&nelem,&ptxyz,&elems);
        if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (elems/ptxyz) failed.\n");
        return 1;
    }
    // created required data structure for mesh 
    int Nredge = 3;
    int *esurp,*esurp_ptr,*esure,*open;
    save_esurp(npoin,nelem,elems,&esurp,&esurp_ptr,Nredge);
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (esurp/esurp_ptr) failed.\n");
        return 1;
    }
    save_esure(nelem,elems,esurp_ptr,esurp,&esure,&open,Nredge);
    if (open == NULL || esure == NULL)
    {
        fprintf(stderr,"Memory allocation (esure/open) failed.\n");
        return 1;
    }
    int opencount = 0;
    for (int i=0;i<nelem;i++){
        if (open[i]==0) opencount++;
    }
    (opencount==0) ? printf("! this is open mesh.\n") : printf("* this is close mesh.\n");
    //calculate the centeroid of each element//
    double *cen;
    CHECK_ERROR(save_centri3(nelem, elems, ptxyz, &cen));
    if (cen == NULL)
    {
        fprintf(stderr, "Memory allocation (cen) failed.\n");
        return 1;
    }
    // calc norm of ele
    double *normele;
    CHECK_ERROR(save_normele(nelem,elems,ptxyz,&normele));
    // flip the normal vector to be outward:
    for (int ele = 0; ele < (3 * nelem); ele++)
        normele[ele] = -1 * normele[ele];
    // read centroid of inlet and outlet hole from inlet.txt and outlet.txt
    FILE *fptr = fopen ("inlet.txt","r");
    if (fptr == NULL){
        fprintf(stderr,"ERROR: can not open inlet.txt file!\n");
        exit(EXIT_FAILURE);
    }
    int buffer = 50, nscan=0;
    char line[buffer];
    int num_inlet=0; 
    double *cen_inlet = (double *)malloc(300*sizeof(double));
    while (fgets(line, buffer, fptr) != NULL) {
        nscan = sscanf(line,"%lf %lf %lf",&cen_inlet[num_inlet*3],&cen_inlet[num_inlet*3+1],&cen_inlet[num_inlet*3+2]);
        if (nscan !=3){
            fprintf(stderr, "ERROR: Invalid file format.\n");
            printf("n=%d string = %s\n",nscan,line);
            break;
        }
        num_inlet++;
    };
    if (num_inlet==0){
        fprintf (stderr, "ERROR: there is no inlet for this case.\n");
        exit(EXIT_FAILURE);
    }
    fclose(fptr);
    #ifdef DEBUG
    for (int i=0;i<num_inlet;i++) printf(" centroid of inlet %i : %lf %lf %lf\n",i,cen_inlet[3*i],cen_inlet[3*i+1],cen_inlet[3*i+2]);
    #endif
    fptr = fopen ("outlet.txt","r");
    if (fptr == NULL){
        fprintf(stderr,"ERROR: can not open outlet.txt file!\n");
        exit(EXIT_FAILURE);
    }
    double *cen_outlet = (double *)malloc(300*sizeof(double));
    int num_outlet = 0;
    while (fgets(line, buffer, fptr) != NULL) {
        nscan = sscanf(line,"%lf %lf %lf",&cen_outlet[num_outlet*3],&cen_outlet[num_outlet*3+1],&cen_outlet[num_outlet*3+2]);
        if (nscan !=3){
            fprintf(stderr, "ERROR: Invalid file format.\n");
            printf("n=%d string = %s\n",nscan,line);
            break;
        }
        num_outlet++;
    };
    if (num_outlet==0){
        fprintf (stderr, "ERROR: there is no outlet for this case.\n");
        exit(EXIT_FAILURE);
    }
    fclose(fptr);
    int num_all_holes = num_inlet + num_outlet;
    double *cen_all_holes = (double *)malloc((size_t)num_all_holes*3*sizeof(double));
    for (int i=0;i<(3*num_inlet);i++) cen_all_holes[i]=cen_inlet[i]; 
    for (int i=0;i<(3*num_outlet);i++) cen_all_holes[3*num_inlet+i]=cen_outlet[i];
    #ifdef DEBUG
    for (int i=0;i<num_outlet;i++) printf("- centroid of outlet %i : %lf %lf %lf\n",i,cen_outlet[3*i],cen_outlet[3*i+1],cen_outlet[3*i+2]);
    #endif
    // mark the clossest element to centroid of hole 
    int *hole_mask = (int *)calloc((size_t)nelem,sizeof(int));
    for (int hole_id=0;hole_id<num_all_holes;hole_id++){
        double min_dist = 100;
        int min_id = -1;
        for (int ele =0;ele<nelem;ele++){
            double ele_cen [3] = {cen[3*ele],cen[3*ele+1],cen[3*ele+2]};
            double dist = SQUARE(ele_cen[0]-cen_all_holes[3*hole_id]); 
            dist += SQUARE(ele_cen[1]-cen_all_holes[3*hole_id+1]); 
            dist += SQUARE(ele_cen[2]-cen_all_holes[3*hole_id+2]); 
            dist = sqrt (dist);
            if (min_dist>dist){
                min_dist = dist;
                min_id = ele;
            }
        } 
        if (min_id == -1){
            fprintf(stderr, "ERROR: No closest element found for hole %d\n", hole_id);
            exit(EXIT_FAILURE);
        }
        //hole_mask [min_id] = hole_id+1;    
        // find all element of each holes 
            int list_capacity = 100; 
            int *list = (int *)calloc((size_t)list_capacity, sizeof(int));
            if (!list) {
                fprintf(stderr, "Memory allocation failed!\n");
                exit(EXIT_FAILURE);
            }
            list [0] = min_id;
            double norm_this_hole[3]={normele[3*min_id],normele[3*min_id+1],normele[3*min_id+2]};
            int list_size = 1;
            do{
                int ele=list[list_size-1];
                double product = norm_this_hole[0]*normele[3*ele]+norm_this_hole[1]*normele[3*ele+1]+norm_this_hole[2]*normele[3*ele+2];
                if (fabs(product)> 0.98){
                    hole_mask [ele] = hole_id+1;
                    list_size--;
                    if(hole_mask[esure [ele*3]]==0) {
                        list_size++;
                        list[list_size-1] = esure [ele*3];
                    }
                    if(hole_mask[esure [ele*3+1]]==0) {
                        list_size++;
                        list[list_size-1] = esure [ele*3+1];
                    }
                    if(hole_mask[esure [ele*3+2]]==0) {
                        list_size++;
                        list[list_size-1] = esure [ele*3+2];
                    }
                }else{
                    list_size--;
                }
                if (list_size >= list_capacity) {
                    list_capacity *= 2;  // Double the capacity
                    list = (int *)realloc(list, (size_t)list_capacity * sizeof(int));
                    if (!list) {
                        fprintf(stderr, "Memory reallocation failed!\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }while(list_size!=0);
            free (list);
    }
    // modified the mask by files 
    int *mask_value = calloc ((size_t)num_all_holes+1,sizeof(int));
    fptr = fopen ("mask_correction.txt","r");
    // Create a hash table
    HashTable table = {0};
    if (fptr != NULL){
        char line_t[100], *str_t, key[50], value[50];
        while ((str_t = fgets(line_t, 100, fptr)))
        {
            if (line_t[0] != '/') {
                char *token = strtok(str_t, " \t\n");
                if (token) {
                    strcpy(key, token);
                    token = strtok(NULL, " \t\n");
                    if (token) {
                        strcpy(value, token);
                        inserthash(&table, key, value);
                    }
                }
            }
        }
        // free file pointer
        fclose(fptr);
        for (int hole_id =0;hole_id<num_all_holes;hole_id++){
            char str[20];  // Make sure the buffer is large enough
            sprintf(str, "%d", hole_id);  // Convert integer to string
            char *hash_result = gethash(&table, str);
            if (hash_result) {
                mask_value[hole_id] = atoi(hash_result);
            } else {
                mask_value[hole_id] = 0;  // Or any default value
            }
        }
        // Free the memory allocated for the hash table
	    freeTable(&table);
    }else{
        fprintf(stderr,"WARNING: save Boundary Condition mask without correction file!\n");
    }
    #ifdef DEBUG
    for(int i=0;i<num_all_holes;i++)
    printf("- value of hole id %d is %d\n",i,mask_value[i]);
    #endif
    for (int ele=0;ele<nelem;ele++){
        if (hole_mask[ele] != 0 ) hole_mask[ele]=mask_value[hole_mask[ele]];
    }
    // creat mesh struct 
    mesh *M1 = (mesh *)malloc(sizeof(mesh));
    if (M1)
    {
        *M1 = (mesh){0}; // Set all integer and pointer fields to 0 or NULL
    }
    if (M1 == NULL)
    {
        fprintf(stderr, "Memory allocation failed for M1 pointer\n");
        exit(EXIT_FAILURE);
    }
    strcpy(M1->type, "tri");
    M1->nrpts = 3;
    M1->elems =elems;M1->npoin =npoin;M1->nelem=nelem;M1->esure=esure;M1->esurp=esurp;M1->esurp_ptr=esurp_ptr;M1->nredge=Nredge;M1->open=open;M1->ptxyz=ptxyz;
    
    // save all data of mesh and data struture in VTK format
    FunctionWithArgs prtelefield[] =
        {
            {"open", 1, nelem, open, SCA_int_VTK},
            {"bc_mask", 1, nelem, hole_mask, SCA_int_VTK},
        };
    size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
    FunctionWithArgs prtpntfield[] = {0};
    size_t countpnt = 0;
    CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
    // allocated memory
    free(elems);
    free(ptxyz);
    free (esure);
    free (esurp);
    free (esurp_ptr);
    free (open);
    free (normele);
    free (M1);
    free (cen_inlet);
    free (cen_outlet);
    free (cen_all_holes);
    free (hole_mask);
    free (cen);
    free (mask_value);
    return 0;
}