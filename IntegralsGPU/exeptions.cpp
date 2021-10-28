#include "exeptions.h"

std::ostream& operator<<(std::ostream& out, const Exeption& exeption)
{
   out << exeption.message;
   return out;
}

ParsingExeption::ParsingExeption() : Exeption("Parsing exeption!") {}

MallocExeption::MallocExeption() : Exeption("Malloc exeption!") {}

CopyExeption::CopyExeption() : Exeption("Copy exeption!") {}
