#pragma once

#include "nlohmann/json.hpp"

#include <string>
#include <fstream>
#include <iomanip>
#include <ios>

using json = nlohmann::json;

namespace tsbp
{

class Serializer
{
public:
    template <class T>
    static void WriteToJson(T& classToSerialize, const std::string& filePath)
    {
        std::ofstream ofs(filePath + ".json");
        {
            if (!ofs.is_open())
            {
                auto ex = std::ios_base::failure("File stream not open.");
                throw std::exception(ex);
            }

            json jsonFile = classToSerialize;
            ofs << std::setw(2) << jsonFile;
        }

        ofs.close();
    }

    template <class T>
    static T ReadFromJson(const std::string& filePath)
    {
        std::ifstream inputStream(filePath);
        nlohmann::json jsonParameters;
        inputStream >> jsonParameters;

        T jsonAsClass(jsonParameters.get<T>());

        return jsonAsClass;
    }
};
    
}