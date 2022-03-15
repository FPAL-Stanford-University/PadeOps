function create_files {
    echo "============================="
    echo "Creating Cdiff_ge input files"
    echo "============================="
    cp input_baseline.dat input_Cge1e-2.dat 
    cp input_baseline.dat input_Cge1e-1.dat 
    cp input_baseline.dat input_Cge1e0.dat 
    cp input_baseline.dat input_Cge1e1.dat 
    cp input_baseline.dat input_Cge1e2.dat 
    cp input_baseline.dat input_Cge1e3.dat 
    
    sed -i 's/out_baseline/Cge\/out_Cge1e-2/g'  input_Cge1e-2.dat
    sed -i 's/out_baseline/Cge\/out_Cge1e-1/g'  input_Cge1e-1.dat
    sed -i 's/out_baseline/Cge\/out_Cge1e0/g'  input_Cge1e0.dat
    sed -i 's/out_baseline/Cge\/out_Cge1e1/g'  input_Cge1e1.dat
    sed -i 's/out_baseline/Cge\/out_Cge1e2/g'  input_Cge1e2.dat
    sed -i 's/out_baseline/Cge\/out_Cge1e3/g'  input_Cge1e3.dat
    
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d-2/g'  input_Cge1e-2.dat
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d-1/g'  input_Cge1e-1.dat
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d0/g'  input_Cge1e0.dat
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d1/g'  input_Cge1e1.dat
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d2/g'  input_Cge1e2.dat
    sed -i 's/Cdiff_g = 0.0d0/Cdiff_g = 1.0d3/g'  input_Cge1e3.dat

    echo "============================="
    echo "Creating Cbeta input files"
    echo "============================="
    cp input_baseline.dat input_Cbeta1e-2.dat 
    cp input_baseline.dat input_Cbeta1e-1.dat 
    cp input_baseline.dat input_Cbeta1e0.dat 
    cp input_baseline.dat input_Cbeta1e1.dat 
    cp input_baseline.dat input_Cbeta1e2.dat 
    cp input_baseline.dat input_Cbeta1e3.dat 
    
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e-2/g'  input_Cbeta1e-2.dat
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e-1/g'  input_Cbeta1e-1.dat
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e0/g'  input_Cbeta1e0.dat
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e1/g'  input_Cbeta1e1.dat
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e2/g'  input_Cbeta1e2.dat
    sed -i 's/out_baseline/Cbeta\/out_Cbeta1e3/g'  input_Cbeta1e3.dat
    
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d-2/g'  input_Cbeta1e-2.dat
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d-1/g'  input_Cbeta1e-1.dat
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d0/g'  input_Cbeta1e0.dat
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d1/g'  input_Cbeta1e1.dat
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d2/g'  input_Cbeta1e2.dat
    sed -i 's/Cbeta = 0.0d0/Cbeta = 1.0d3/g'  input_Cbeta1e3.dat

    echo "============================="
    echo "Setting LAD variations"
    echo "============================="
    sed -i 's/g_LAD_id = 0/g_LAD_id = 0/g'    input_C*.dat
    sed -i 's/gp_LAD_id = 0/gp_LAD_id = 0/g'  input_C*.dat
    sed -i 's/gt_LAD_id = 0/gt_LAD_id = 0/g'  input_C*.dat
    sed -i 's/beta_LAD_id = 0/beta_LAD_id = 1/g'  input_C*.dat

    echo "============================="
    echo "============================="
}
