#include "exeptions.hpp"

std::ostream& operator<<(std::ostream& out, const Exeption& exeption)
{
   out << exeption.message;
   return out;
}

Exeption::Exeption() : message("Error!") {}
Exeption::Exeption(std::string message) : message(message) {}

ParsingExeption::ParsingExeption() : Exeption("Parsing exeption!") {}

MallocExeption::MallocExeption() : Exeption("Malloc exeption!") {}

CopyExeption::CopyExeption() : Exeption("Copy exeption!") {}
