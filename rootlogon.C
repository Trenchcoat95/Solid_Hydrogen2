{
    // Print confirmation (optional)
    std::cout << "Loading SAND libraries..." << std::endl;

    // Path to your compiled libraries (adjust as needed)
    const char* lib_path = "./";
    

    // Load ROOT dependencies first (if needed)
    gSystem->Load("libPhysics");
    gSystem->Load("libTree");

    // Load all .so libraries (explicitly list them)
    const std::vector<TString> libraries = {
        "../sandreco/lib/libStruct.so",
        "./My_libMyDict.so"
    };

    for (const auto& lib : libraries) {
        if (gSystem->Load(lib_path + lib) < 0) {
            std::cerr << "Error loading " << lib << std::endl;
        }
    }

    if(gSystem->Load("/opt/exp_software/neutrino/al9/EDEPSIM/edep-sim/lib/libedepsim_io.so")<0){
       std::cerr << "Error loading EDEPSIM" << std::endl;

      }

    

    // Print all loaded libraries
    std::cout << "\n=== Loaded Libraries ===\n";
    auto libs = gSystem->GetLibraries();
    std::cout<< libs<<"\n";
    std::cout << "======================\n";

}

