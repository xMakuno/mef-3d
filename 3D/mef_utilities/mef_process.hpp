float calculate_local_area(float x1, float y1, float x2, float y2, float x3, float y3){
    float A = abs((x1*y2 + x2*y3 + x3*y1) - (x1*y3 + x2*y1 + x3*y2))/2;
    return ((A==0)?0.000001:A);
}

//TODO: add calculate_local_volume function to calculate an elements volume
//DONE
float calculate_local_volume(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4){
    // This is a 4x4 matrix determinant
    /*
        x1 y1 z1 1    
        x2 y2 z2 1    
        x3 y3 z3 1    
        x4 y4 z4 1

        A
        y2 z2 1    
        y3 z3 1    
        y4 z4 1

        B
        x2 z2 1
        x3 z3 1
        x4 z4 1

        C
        x2 y2 1
        x3 y3 1
        x4 y4 1

        D
        x2 y2 z2    
        x3 y3 z3    
        x4 y4 z4
    */
   /* float detA = y2*z3 + y3*z4 + y4*z2 - y4*z3 - y2*z4 - y3*z2;
   float detB = x2*z3 + x3*z4 + x4*z2 - x4*z3 - x2*z4 - x3*z2;
   float detC = x2*y3 + x3*y4 + x4*y2 - x4*y3 - x2*y4 - x3*y2; 
   float detD = x2*y3*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2 - x2*y4*z3 - x3*y2*z4;
   */
   float detA = y2*(z3-z4) - z2*(y3-y4) + ((y3*z4)-(y4*z3));
   float detB = x2*(z3-z4) - z2*(x3-x4) + ((x3*z4)-(x4*z3));
   float detC = x2*(y3-y4) - y2*(x3-x4) + ((x3*y4)-(x4*y3));
   float detD = x2*((y3*z4)-(y4*z3)) - y2*((x3*z4)*(x4*z3)) + z2*((x3*y4)-(x4*y3));
   float V = abs(x1*detA - y1*detB* + z1*detC - detD)/6;
   
   return ((V==0)?0.000001:V);
}

//TODO: update local jacobian function to include the Z axis
// DONE
float calculate_local_jacobian(float x1, float y1, float z1, float x2,  float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4){
    //TODO: Change jacobian formula to insanely big formula
    // DONE
    /*
        x21 x31 x41
        y21 y31 y41
        z21 z31 z41

        A
        y31 y41
        z31 z41
        B
        y21 y41
        z21 z41
        C
        y21 y31
        z21 z31
    */
    float detA = (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1);
    float detB = (y2-y1)*(z4-z1) - (y4-y1)*(z2-z1);
    float detC = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
    float J = (x2-x1)*detA + (x1-x3)*detB + (x4-x1)*detC;
    // float J = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
    return ((J==0)?0.000001:J);
}

//TODO: refactor calculate B to meet the following:
/*
    -1 1 0 0
    -1 0 1 0
    -1 0 0 1
 */
// DONE
void calculate_B(Matrix* B){
    B->set(-1,0,0);  B->set(1,0,1);  B->set(0,0,2);  B->set(0,0,3);
    B->set(-1,1,0);  B->set(0,1,1);  B->set(1,1,2);  B->set(0,1,3);
    B->set(-1,2,0);  B->set(0,2,1);  B->set(0,2,2);  B->set(1,2,3);
}

//TODO: change function to meet the dimensions of new A matrix in 3d
// DONE
void calculate_local_A(Matrix* A, float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4){
    /* A->set(y3-y1, 0, 0);   A->set(x1-x3, 0, 1);  
    A->set(y1-y2, 1, 0);   A->set(x2-x1, 1, 1); */
    float a = (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1);
    float b = (x4-x1)*(z3-z1) - (x3-x1)*(z4-z1);
    float c = (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1);
    float d = (y4-y1)*(z2-z1) - (y2-y1)*(z4-z1);
    float e = (x2-x1)*(z4-z1) - (x4-x1)*(z2-z1);
    float f = (x4-x1)*(y2-y1) - (x2-x1)*(y4-y1);
    float g = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1); 
    float h = (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1);
    float i = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    A->set(a,0,0);A->set(d,0,1);A->set(g,0,2);
    A->set(b,1,0);A->set(e,1,1);A->set(h,1,2);
    A->set(c,2,0);A->set(f,2,1);A->set(i,2,2);
}

//TODO: change function to implement all changes to new K matrices
void create_local_K(Matrix* K, int element_id, Mesh* M){
    //TODO: change dimensions to 4 by 4
    // DONE
    K->set_size(4,4);

    float k = M->get_problem_data(THERMAL_CONDUCTIVITY);
    //TODO: add new coordinates for K matrix
    // DONE
    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();
    //TODO: update to implemente volume function
    // DONE
    // float Area = calculate_local_area(x1, y1, x2, y2, x3, y3);
    float Volume = calculate_local_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    cout << "\tVolume: " << Volume << endl;
    float J = calculate_local_jacobian(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    
    //TODO: change matrices declaration to meet new dimensions
    // DONE
    Matrix B(3,4), A(3,3);
    calculate_B(&B);
    //TODO: update number of arguments with the new coordinates
    calculate_local_A(&A, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    //B.show(); A.show();

    //TODO: change dimensions in these matrices to meet the new dimensions
    //DONE
    Matrix Bt(4,3), At(3,3);
    transpose(&B,3,4,&Bt);
    transpose(&A,3,3,&At);
    //Bt.show(); At.show();

    Matrix res1, res2, res3;
    product_matrix_by_matrix(&A,&B,&res1);
    product_matrix_by_matrix(&At,&res1,&res2);
    product_matrix_by_matrix(&Bt,&res2,&res3);
    // TODO: Update formula to pass in the Volume & update dimensions
    // DONE
    product_scalar_by_matrix(k*Volume/(J*J),&res3,4,4,K);

    //cout << "\t\tLocal matrix created for Element " << element_id+1 << ": "; K->show(); cout << "\n";
}

void create_local_b(Vector* b, int element_id, Mesh* M){
    //TODO: update local B size to 4
    // DONE
    b->set_size(4);

    float Q = M->get_problem_data(HEAT_SOURCE);
    //TODO: add new coorddinates x4 and y4 and all z coordinates
    // DONE
    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();
    //TODO: update the parameters that are passed in the function
    // DONE
    float J = calculate_local_jacobian(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

    // TODO: update b formula to Q*J / 24 and add a 4th row in local b matrix
    // DONE
    b->set(Q*J/24.0f, 0);
    b->set(Q*J/24.0f, 1);
    b->set(Q*J/24.0f, 2);
    b->set(Q*J/24.0f, 3);

    //cout << "\t\tLocal vector created for Element " << element_id+1 << ": "; b->show(); cout << "\n";
}

void create_local_systems(Matrix* Ks, Vector* bs, int num_elements, Mesh* M){
    for(int e = 0; e < num_elements; e++){
        cout << "\tCreating local system for Element " << e+1 << "...\n\n";
        create_local_K(&Ks[e],e,M);
        create_local_b(&bs[e],e,M);
    }
}

//TODO: add new row and column to meet new dimensions of 4
//TODO: add index4 to parameters
// DONE
void assembly_K(Matrix* K, Matrix* local_K, int index1, int index2, int index3, int index4){
    K->add(local_K->get(0,0),index1,index1);    K->add(local_K->get(0,1),index1,index2);    K->add(local_K->get(0,2),index1,index3);    K->add(local_K->get(0,3),index1,index4);
    K->add(local_K->get(1,0),index2,index1);    K->add(local_K->get(1,1),index2,index2);    K->add(local_K->get(1,2),index2,index3);    K->add(local_K->get(1,3),index2,index4);
    K->add(local_K->get(2,0),index3,index1);    K->add(local_K->get(2,1),index3,index2);    K->add(local_K->get(2,2),index3,index3);    K->add(local_K->get(2,3),index3,index4);
    K->add(local_K->get(3,0),index4,index1);    K->add(local_K->get(3,1),index4,index2);    K->add(local_K->get(3,2),index4,index3);    K->add(local_K->get(3,3),index4,index4);
}

//TODO: add 4th row to calculation
// DONE
void assembly_b(Vector* b, Vector* local_b, int index1, int index2, int index3, int index4){
    b->add(local_b->get(0),index1);
    b->add(local_b->get(1),index2);
    b->add(local_b->get(2),index3);
    b->add(local_b->get(3),index4);
}


void assembly(Matrix* K, Vector* b, Matrix* Ks, Vector* bs, int num_elements, Mesh* M){
    K->init();
    b->init();
    //K->show(); b->show();
    //TODO: add 4th node to assembly
    for(int e = 0; e < num_elements; e++){
        cout << "\tAssembling for Element " << e+1 << "...\n\n";
        int index1 = M->get_element(e)->get_node1()->get_ID() - 1;
        int index2 = M->get_element(e)->get_node2()->get_ID() - 1;
        int index3 = M->get_element(e)->get_node3()->get_ID() - 1;
        int index4 = M->get_element(e)->get_node4()->get_ID() - 1;
        // TODO: pass index4 to functions
        assembly_K(K, &Ks[e], index1, index2, index3, index4);
        assembly_b(b, &bs[e], index1, index2, index3, index4);
        //cout << "\t\t"; K->show(); cout << "\t\t"; b->show(); cout << "\n";
    }
}

void apply_neumann_boundary_conditions(Vector* b, Mesh* M){
    int num_conditions = M->get_quantity(NUM_NEUMANN);

    for(int c = 0; c < num_conditions; c++){
        Condition* cond = M->get_neumann_condition(c);
        
        int index = cond->get_node()->get_ID() - 1;
        b->add(cond->get_value(), index);
    }
    //cout << "\t\t"; b->show(); cout << "\n";
}

void add_column_to_RHS(Matrix* K, Vector* b, int col, float T_bar){
    for(int r = 0; r < K->get_nrows(); r++)
        b->add(-T_bar*K->get(r,col),r);
}

void apply_dirichlet_boundary_conditions(Matrix* K, Vector* b, Mesh* M){
    int num_conditions = M->get_quantity(NUM_DIRICHLET);
    int previous_removed = 0;

    for(int c = 0; c < num_conditions; c++){
        Condition* cond = M->get_dirichlet_condition(c);
        
        int index = cond->get_node()->get_ID() - 1 - previous_removed;
        float cond_value = cond->get_value();

        //K->show();
        K->remove_row(index);
        //K->show();
        //b->show();
        b->remove_row(index);
        //b->show();

        add_column_to_RHS(K, b, index, cond_value);
        //b->show();

        K->remove_column(index);
        //K->show();

        previous_removed++;
    }
}

void solve_system(Matrix* K, Vector* b, Vector* T, int mode){
    int n = K->get_nrows();
    
    Matrix Kinv(n,n);

    cout << "\tCalculating inverse of global matrix K...\n\n";
    if(mode == 1) calculate_inverse(K, n, &Kinv); //1
    else calculate_inverse_Cholesky(K, n, &Kinv); //2

    cout << "\tPerforming final calculation...\n\n";
    product_matrix_by_vector(&Kinv, b, n, n, T);
}

void merge_results_with_dirichlet(Vector* T, Vector* Tf, int n, Mesh* M){
    int num_dirichlet = M->get_quantity(NUM_DIRICHLET);

    int cont_dirichlet = 0;
    int cont_T = 0;
    for(int i = 0; i < n; i++){
        if(M->does_node_have_dirichlet_condition(i+1)){
            Condition* cond = M->get_dirichlet_condition(cont_dirichlet);
            cont_dirichlet++;
        
            float cond_value = cond->get_value();

            Tf->set(cond_value,i);
        }else{
            Tf->set(T->get(cont_T),i);
            cont_T++;
        }
    }
}
