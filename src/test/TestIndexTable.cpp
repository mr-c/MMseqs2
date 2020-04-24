//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//

#include <cstdio>
#include <iostream>
#include "SubstitutionMatrix.h"
#include "IndexTable.h"
#include "IndexBuilder.h"
#include "Parameters.h"

const char* binary_name = "test_indextable";

int main (int, const char**) {
    Parameters &par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 8.0, -0.2f);
    DBReader<unsigned int> dbr(
                               "../example-data/DB", "../example-data/DB.index", 1, 1);
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Sequence *s = new Sequence(32000, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 6, true, false);
    IndexTable t(subMat.alphabetSize, 6, false);
    IndexBuilder::fillDatabase(&t, NULL, NULL, subMat, s, &dbr, 0, dbr.getSize(), 0, 1, 1);
    t.printStatistics(subMat.num2aa);

    delete s;
    dbr.close();

    return 0;
}

