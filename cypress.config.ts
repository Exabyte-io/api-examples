import { defineConfig } from 'cypress'

export default defineConfig({
  e2e: {
    baseUrl: 'https://jupyterlite.readthedocs.io/en/latest/_static/lab',
    supportFile: 'cypress/support/e2e.ts',
    specPattern: 'cypress/e2e/**/*.cy.{js,jsx,ts,tsx}',
    video: false,
    screenshotOnRunFailure: false,
    defaultCommandTimeout: 20000,
  },
}) 