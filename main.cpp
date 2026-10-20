#include <iostream>

#include "Config.h"
#include "Applications/Applications.h"

int main(int argc, char* argv[]){

    if(argc < 2){
        std::cout << "To start NASTS, USE THIS COMMAND: NASTS *input_file_name.json*\n";
        return 1;
    }

    //Read and parse the input json file provided by the user and store it in RAM
    Config::access().read(argv[1]);

    //See the requested solver type
    std::string appType = Config::access().data()["Application"];

    //run the requested solver
    auto solver = Applications::make(appType);
    solver->run();

    return 0;
}