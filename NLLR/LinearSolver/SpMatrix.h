// SpMatrix.h: interface for the CSpMatrix class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _LINEAR_SOLVER_SP_MATRIX_H_INCLUDED_
#define _LINEAR_SOLVER_SP_MATRIX_H_INCLUDED_

class CSpMatrix	// Intermediate Sparse Matrix
{
	// Data Conversion
	friend class CTAUCS;
	friend class CSuperLU;

public:
	CSpMatrix(const int rowNum = 0, const int colNum = 0);	// Default Constructor
	CSpMatrix(const CSpMatrix &src);						// Copy Constructor
	virtual CSpMatrix& operator=(const CSpMatrix &src);		// Assignment Operator
	virtual ~CSpMatrix();									// Destructor
	int RowNum() const {	return m_rowNum;	};	// Number of Rows
	int ColNum() const {	return m_colNum;	};	// Number of Columns
	int NonZeroElementNum() const {	return m_nnzNum;	}; // Number of Non-zero Entries // Added by Yunbo 2012-12-04 
	void SetRowCol(bool flag) {	  m_bRow = flag;   };
	bool Size(const int rowNum, const int colNum);	// Re-size and clear all data
	bool Set(const int iRow, const int iCol, const double value, const bool bAdd = false);	// Set or add value to entry
	bool Get(const int iRow, const int iCol, double &value) const;							// Get value of entry
	static bool Transpose(const CSpMatrix &in, CSpMatrix &out);	// Transpose
	static bool MulATA(CSpMatrix &matA);				// Multiplication: A'A
	static bool MulATA(CSpMatrix &matA, double *&b);	// Multiplication: Ax = b => A'Ax = A'b

	static bool MulATA(CSpMatrix &matA, CSpMatrix& matATA);
	static bool MulATB(CSpMatrix &matA, double* b, double* ATb);	// Multiplication: ATb = A'b

protected:
	typedef struct	// Row or Column
	{
		int m_num, *m_ids;	// No. and Indices of Columns or Rows with Non-zero Enties
		double *m_vals;		// Values of Non-zero Entries
	} RowCol;

	// Functions
	bool Switch();	// Switch to row- or column-wise

	// Variables
	bool m_bRow;						// Whether storage is row- or column-wise
	int m_rowNum, m_colNum, m_nnzNum;	// No. of Rows, Columns and Non-zero Entries
	RowCol *m_rowCols;					// Array of Rows or Columns
};

#endif	// _LINEAR_SOLVER_SP_MATRIX_H_INCLUDED_
