#pragma once
#include <string>
#include <ostream>

class NotImplementedExeption{};
class RangeExeption{};

class Exeption
{
public:
   std::string message;
   Exeption();
   Exeption(std::string message);
};

class ParsingExeption : public Exeption
{
public:
   ParsingExeption();
};

class MallocExeption : public Exeption
{
public:
   MallocExeption();
};

class CopyExeption : public Exeption
{
public:
   CopyExeption();
};

std::ostream& operator <<(std::ostream& out, const Exeption& exeption);
