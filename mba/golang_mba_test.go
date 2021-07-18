package mba

import (
	"math"
	"testing"
)

func TestMba2Trivial(t *testing.T) {
	m, err := NewMba2([]float64{0, 0}, []float64{10, 10}, []int64{3, 3},
		[]float64{5, 6},
		[]float64{7, 8},
		[]float64{3, 4},
		nil)
	if err != nil {
		t.Error(err)
	}
	defer m.Free()

	// fmt.Println(m.String())
}

func TestMba2KnownValues(t *testing.T) {
	// based on the python parameters in the README

	interp, err := NewMba2([]float64{-0.1, -0.1},
		[]float64{1.1, 1.1},
		[]int64{3, 3},
		[]float64{0.0, 0.0, 1.0, 1.0, 0.4, 0.6},
		[]float64{0.0, 1.0, 0.0, 1.0, 0.4, 0.6},
		[]float64{0.2, 0.0, 0.0, -0.2, -1.0, 1.0},
		nil)
	if err != nil {
		t.Fatal(err)
	}
	defer interp.Free()

	x := []float64{0.1, 0.7, 0.4, 0.4, 0.5}
	y := []float64{0.7, 0.7, 0.7, 0.6, 0.6}

	// test Apply
	val, err := interp.Apply(x, y)
	if err != nil {
		t.Error(err)
	}
	expected := []float64{-0.20030558288394174, 0.8637352407147169, 0.3290533037278714, 5.305814381600937e-16, 0.5470065204308358}

	if len(expected) != len(val) {
		t.Errorf("Return vector from Apply has unexpected length %d, expected %d", len(val), len(expected))
	}

	for i := 0; i < len(val); i++ {
		if math.Abs(val[i]-expected[i]) > 1e-8 {
			t.Errorf("Exceeded tolerance in Apply: expected value %f got %f", expected[i], val[i])
		}
	}

	// test ApplySingle
	for i := 0; i < len(x); i++ {
		val := interp.ApplySingle(x[i], y[i])
		if math.Abs(val-expected[i]) > 1e-8 {
			t.Errorf("Exceeded tolerance in ApplySingle: expected value %f got %f", expected[i], val)
		}
	}
}

func TestMba2(t *testing.T) {

	// From MBA tests's "small test function"
	coo_x := []float64{}
	coo_y := []float64{}
	val := []float64{}

	for y := 0; y < 10; y++ {
		for x := 0; x < 10; x++ {
			if (x%2 == 0) && (y%2 == 0) {
				x_coord := float64(x) / 10.
				y_coord := float64(y) / 10.
				current_value := math.Sin(float64(x)*8.) * math.Sin(float64(y)*8.)
				coo_x = append(coo_x, x_coord)
				coo_y = append(coo_y, y_coord)
				val = append(val, current_value)
			}
		}
	}

	lo := []float64{-0.1, -0.1}
	hi := []float64{1.1, 1.1}
	grid := []int64{2, 2}
	phi, err := NewMba2(lo, hi, grid, coo_x, coo_y, val, &Config{MaxLevels: 8, Tol: 1e-8, MinFill: 0.5})
	if err != nil {
		t.Fatal(err)
	}
	defer phi.Free()

	// apply to the coordinates
	computed, err := phi.Apply(coo_x, coo_y)
	if err != nil {
		t.Errorf("Error in bulk apply: %s", err)
	}
	if len(computed) != len(val) {
		t.Errorf("Return vector from Apply has unexpected length %d, expected %d", len(computed), len(val))
	}
	for i := 0; i < len(val); i++ {
		if math.Abs(val[i]-computed[i]) > 1e-8 {
			t.Errorf("Exceeded tolerance: coordinate (%f, %f) expected value %f got %f", coo_x[i], coo_y[i], val[i], computed[i])
		}
	}

	// check like in the test
	for y := 0; y < 10; y++ {
		for x := 0; x < 10; x++ {
			x_coord := float64(x) / 10.
			y_coord := float64(y) / 10.
			true_value := math.Sin(float64(x)*8.) * math.Sin(float64(y)*8.)
			inter_value, err := phi.Apply([]float64{x_coord}, []float64{y_coord})
			if err != nil {
				t.Errorf("Error in Apply: %s", err)
			}

			if math.Abs(true_value-inter_value[0]) >= 1.5 {
				t.Error("unexpected mismatch")
			}

			// also test ApplySingle
			inter_value_single := phi.ApplySingle(x_coord, y_coord)
			if math.Abs(true_value-inter_value_single) >= 1.5 {
				t.Error("unexpected mismatch")
			}

		}

	}
}
