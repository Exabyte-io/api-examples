import { NotebookTester } from '../../support/notebook-tester'

describe('Materials Designer Specific Examples Notebooks', () => {
  const tester = new NotebookTester()

  beforeEach(() => {
    cy.viewport(1920, 1080)
  })

  const notebookTests = [
    '/other/materials_designer/specific_examples/defect_planar_grain_boundary_2d_boron_nitride.ipynb'
    // Add other notebook paths as needed
  ]

  notebookTests.forEach((notebookPath) => {
    it(`should run ${notebookPath} successfully`, () => {
      // Visit the notebook URL
      tester.visit(notebookPath)

      // Wait for Pyodide to be ready
      tester.waitForPyodide()

      // Run all cells and verify no errors
      tester.runAllCells()

      // Check that no cells have errors
      tester.checkNoErrors()
    })
  })
}) 