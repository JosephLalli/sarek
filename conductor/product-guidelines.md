# Product Guidelines

## Prose Style
*   **Tone:** Formal, Technical, and Precise.
*   **Voice:** Active voice is preferred for instructions ("Run the pipeline..."), but passive voice is acceptable for describing system states ("The file is generated...").
*   **Detail:** Prioritize technical accuracy and completeness. Documentation should assume a competent technical audience (bioinformaticians/developers).
*   **Commits:** Follow the Conventional Commits specification (e.g., `feat:`, `fix:`, `refactor:`) with clear, imperative descriptions.

## Visual Identity
*   **Branding:** Strictly adhere to the standard `nf-core` visual identity.
*   **Logos/Assets:** Use existing nf-core/sarek assets without modification.
*   **Reports:** MultiQC and other reports should maintain the default nf-core styling to ensure familiarity for existing users.

## Coding Standards & "Rules of the Road"
*   **Source of Truth:** Strict adherence to the standards defined in `conductor/code_styleguides/nextflow.md` and `conductor/code_styleguides/nf-core.md`.
*   **Strict Syntax Enforcement:**
    *   **Explicit Closures:** absolutely NO implicit `it` variables in closures; all parameters must be named (e.g., `map { file -> ... }`).
    *   **Channel Syntax:** Use `channel.` (lowercase) exclusively; `Channel.` (capitalized) is forbidden.
    *   **Module Structure:** Local modules must use `task.ext.args` and `task.ext.prefix` for argument handling and output naming.
*   **Documentation:**
    *   **Meta Files:** Every local module MUST have a corresponding `meta.yml` defining inputs, outputs, and authors.
    *   **Comments:** Code should be self-documenting where possible, but complex logic (especially channel manipulations) requires explanatory comments.
