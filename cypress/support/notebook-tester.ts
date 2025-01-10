export class NotebookTester {
  visit(notebookPath: string) {
    const cleanPath = notebookPath.replace(/^\//, '')
    cy.visit(`/index.html?path=${cleanPath}`)
  }

  waitForPyodide() {
    // Wait for Pyodide to be initialized
    cy.window().should('have.property', 'pyodide')
    // Additional wait to ensure kernel is fully ready
    cy.wait(2000)
  }

  runAllCells() {
    // Wait for JupyterLab to be ready
    cy.get('#jupyterlab-splash', { timeout: 20000 }).should('not.exist')
    
    // Click the "Run All Cells" button
    cy.get('[data-command="runmenu:run-all"]').click()
  }

  checkNoErrors() {
    // Check for error outputs in cells
    cy.get('.jp-OutputArea-output.jp-OutputArea-executeResult')
      .should('not.contain', 'Error')
    cy.get('.jp-OutputArea-output.jp-OutputArea-stderrOrError')
      .should('not.exist')
  }
} 