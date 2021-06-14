package mba

// #cgo CXXFLAGS: -std=c++11
// #include "golang_mba.h"
import "C"
import (
	"errors"
	"unsafe"
)

var (
	ErrLoDim                  = errors.New("'lo' parameter length should be the same as the MBA dimension")
	ErrHiDim                  = errors.New("'hi' parameter length should be the same as the MBA dimension")
	ErrGridDim                = errors.New("'grid' parameter length should be the same as the MBA dimension")
	ErrCooValueLengthMismatch = errors.New("coordinates and value array lengths should all be equal")
	ErrZeroLength             = errors.New("coordinates and value need to have non-zero length")
)

// Mba2 performs two-dimensional Multilevel B-Spline fitting and evaluation.
type Mba2 struct {
	ptr unsafe.Pointer
}

// Config controls the leveling
type Config struct {
	MaxLevels int
	Tol       float64
	MinFill   float64
}

// NewMba2 creates a new two-dimensional Multilevel B-Spline
// This struct wraps a C++ pointer: it is the caller's responsibility to call Free() to free up memory once no longer used.
//  Inputs:
//   hi, lo, and grid must be of length 2
//   x, y, and val must be of the same non-zero length
//   if config is nil, default values MaxLevels: 8, Tol: 1e-8, MinFill: 0.5 are filled.
//   if config is non-nil:
//      - if MaxLevels==0, the default MaxLevels is set (default value is 8)
//      - if Tol==0, the default Tol is set (default value 1e-8)
//      - since MinFill==0 is a valid value, it is not updated to the default of 0.5
func NewMba2(lo, hi []float64, grid []int64, x, y, val []float64, config *Config) (Mba2, error) {
	// check input variables
	if len(lo) != 2 {
		return Mba2{}, ErrLoDim
	}
	if len(hi) != 2 {
		return Mba2{}, ErrHiDim
	}
	if len(grid) != 2 {
		return Mba2{}, ErrGridDim
	}
	n := len(x)
	if (len(y) != n) || (len(val) != n) {
		return Mba2{}, ErrCooValueLengthMismatch
	}
	if n == 0 {
		return Mba2{}, ErrZeroLength
	}

	// Set config defaults
	if config == nil {
		config = new(Config)
		config.MinFill = 0.5
	}
	if config.MaxLevels == 0 {
		config.MaxLevels = 8
	}
	if config.Tol == 0 {
		config.Tol = 1e-8
	}
	// MinFill = 0 is a valid value, and while the Python library uses 0.5, we leave the user's value,
	// and only set to 0.5 when passed a nil Config

	m := Mba2{}
	m.ptr = C.Mba_2_New((*C.double)(&lo[0]),
		(*C.double)(&hi[0]),
		(*C.longlong)(&grid[0]),
		(C.int)(n),
		(*C.double)(&x[0]),
		(*C.double)(&y[0]),
		(*C.double)(&val[0]),
		(C.int)(config.MaxLevels),
		(C.double)(config.Tol),
		(C.double)(config.MinFill),
	)
	return m, nil
}

// Free frees the C++ class underlying the given Mba2
func (m Mba2) Free() {
	C.Mba_2_Destroy(m.ptr)
}

// String returns the encoding of the Mba2 into a stream, which contains information about levels and their sizes
func (m Mba2) String() string {
	p := C.Mba_2_String(m.ptr)
	s := C.GoString(p)
	C.free(unsafe.Pointer(p))
	return s
}

// Apply uses the Multilevel B-spline to estimate function values at the given coordinates
// the length of x and y must be equal
func (m Mba2) Apply(x, y []float64) ([]float64, error) {
	n := len(x)
	if len(y) != n {
		return nil, ErrCooValueLengthMismatch
	}

	val := make([]float64, n)
	C.Mba_2_Apply(m.ptr,
		(C.int)(n),
		(*C.double)(&x[0]),
		(*C.double)(&y[0]),
		(*C.double)(&val[0]),
	)
	return val, nil

}
