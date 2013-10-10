int main () {}

double Abs(double Nbr)
{
//	return (Nbr >= 0) ? Nbr : -Nbr;
	if( Nbr >= 0 )
		return Nbr;
	else
		return -Nbr;
}

void setResidual(int *solvAtoms, int rMax, int j, int tot, int *residuals, int *vRes) {
  int i;

  for (i = 1; i < tot; i++) {
    residuals[i] = rMax;
    vRes[i] = i;
    if (! solvAtoms[i] || i == j) vRes[i] = 0;
  }
}

int getResidual(int arrayIndex, int **distArray, int *residuals, int *validRes, int tot, int boxLen) {
  int i, j, atomC, tmp[(tot + 1)], count;
  double baseVal, offset;

  baseVal = distArray[arrayIndex][0];
  count = 0;

  for (i = 1; i <= tot; i++) {
    tmp[i] = validRes[i];
    validRes[i] = 0;
  }

  for (j = (arrayIndex +1); j < tot; j++) {
    atomC = distArray[j][1];
    if (! tmp[atomC]) continue;
    offset = Abs(distArray[j][0] - baseVal);
    while (offset > boxLen) offset -= boxLen;
    residuals[atomC] -= (offset*offset);
    if (residuals[atomC] < 0) break;
    validRes[atomC] = atomC;
    count++;
  }

  for (j = (arrayIndex - 1); j > -1 ; j--) {
    atomC = distArray[j][1];
    if (! tmp[atomC]) continue;
    offset = Abs(distArray[j][0] - baseVal);
    while (offset > boxLen) offset -= boxLen;
    residuals[atomC] -= (offset*offset);
    if (residuals[atomC] < 0) break;
    validRes[atomC] = atomC;
    count++;
  }

  return count;
}
