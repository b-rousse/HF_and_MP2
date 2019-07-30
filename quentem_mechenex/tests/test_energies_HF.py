"""
Unit and regression test for the cookie_cutter_example package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
import numpy as np

@pytest.fixture()
def water_molecule():
    name = 'water'
    symbols = ['H', 'O', 'H']
    coordinates = np.array([[2, 0, 0], [0, 0, 0], [-2, 0, 0]])

    water = cookie_cutter_example.Molecule(name, symbols, coordinates)

    return water

def test_molecule_set_coordinates(water_molecule):
    '''Test that bond list is rebuilt when we reset coordinates.'''

    num_bonds = len(water_molecule.bonds)

    assert num_bonds == 2 #makes sure the number of bonds calculated on the initial structure is correct, ie numbonds method is working as expected.

    new_coordinates = np.array([[5, 0, 0], [0, 0, 0], [-2, 0, 0]])
    water_molecule.coordinates = new_coordinates
    new_num_bonds = len(water_molecule.bonds)

    assert new_num_bonds == 1 #tests that the number of bonds is being recalculated
    assert np.array_equal(new_coordinates, water_molecule.coordinates) #tests that the coordinates indeed are being updated

def test_cookie_cutter_example_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "cookie_cutter_example" in sys.modules

def test_calculate_distance():
    '''test the calculate_distance function'''

    r1 = np.array([0,0,-1])
    r2 = np.array([0,1,0])

    expected_distance = np.sqrt(2.)

    calculated_distance = cookie_cutter_example.calculate_distance(r1, r2)

    assert expected_distance == calculated_distance


def test_calculate_angle_90():
    '''test the calculate_angle function'''

    r1 = np.array([1,0,0])
    r2 = np.array([0,0,0])
    r3 = np.array([0,1,0])
    
    calculated_angle = cookie_cutter_example.calculate_angle(r1, r2, r3, degrees = True)

    expected_angle = 90

    assert expected_angle == calculated_angle

def test_calculate_angle_60():
    '''tests another value of the calculate_angle function'''

    r1 = np.array([0,0,-1])
    r2 = np.array([0,1,0])
    r3 = np.array([1,0,0])

    calculated_angle = cookie_cutter_example.calculate_angle(r1, r2, r3, degrees = True)

    expected_angle = 60

    assert np.isclose(expected_angle, calculated_angle, rtol = 1e-9)

#the decorator below helps us automate the generation of tests. The parametrize decorator loops through the test it decorates and runs tests with the information we've hard-coded into the decorator.
@pytest.mark.parametrize("p1, p2, p3, expected_angle", [#test cases are all gonna follow this test template and have these four variables
    (np.array([1, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 0]), 90),#here are the four things for test 1
    (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60),#here are the four things for test 2.
    #so we've written the decorator.
])
def test_calculate_angle(p1, p2, p3, expected_angle):#this is the first function below the decorator, so it is decorated. None that follow will be decorated.

    calculated_angle = cookie_cutter_example.calculate_angle(p1, p2, p3, degrees=True)

    assert np.isclose(expected_angle, calculated_angle, rtol = 1e-9)


def test_create_failure():
    name = 25
    symbols = ['H', 'O', 'H']
    coordinates = np.zeros([3,3])

    with pytest.raises(TypeError): #raises is a way to tell pytest we expect the code to raise an error of type in parentheses
        water = cookie_cutter_example.Molecule(name,symbols,coordinates)
