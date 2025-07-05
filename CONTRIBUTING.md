# Contributing to TMM

Thank you for your interest in contributing to the TMM (Transfer Matrix Method) package! This guide outlines how to contribute effectively to the project.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Contribution Guidelines](#contribution-guidelines)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Community Guidelines](#community-guidelines)

## Getting Started

### Types of Contributions

We welcome various types of contributions:

- **Bug reports**: Help us identify and fix issues
- **Feature requests**: Suggest new features or improvements
- **Code contributions**: Implement new features or fix bugs
- **Documentation**: Improve documentation, tutorials, or examples
- **Testing**: Add or improve test coverage
- **Performance improvements**: Optimize existing code

### Before You Start

1. Check existing [issues](https://github.com/dingxu/tmm-optics/issues) and [pull requests](https://github.com/dingxu/tmm-optics/pulls)
2. For major changes, open an issue first to discuss your approach
3. Read through the codebase to understand the architecture
4. Review this contributing guide and our [Code of Conduct](CODE_OF_CONDUCT.md)

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- A text editor or IDE (we recommend VSCode, PyCharm, or similar)

### Setting Up the Development Environment

1. **Fork and clone the repository**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/tmm-optics.git
   cd tmm-optics
   ```

2. **Create a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install development dependencies**:
   ```bash
   pip install -e ".[dev]"
   ```

4. **Verify the installation**:
   ```bash
   python -c "import tmm; print(tmm.__version__)"
   pytest tests/ -v
   ```

### Development Tools

We use several tools to maintain code quality:

- **Black**: Code formatting
- **Flake8**: Code linting
- **MyPy**: Type checking
- **Pytest**: Testing framework
- **Pre-commit**: Git hooks for automated checks

Install pre-commit hooks:
```bash
pre-commit install
```

## Contribution Guidelines

### Branching Strategy

- `main`: Stable release branch
- `develop`: Integration branch for new features
- `feature/feature-name`: Feature development branches
- `bugfix/bug-description`: Bug fix branches
- `hotfix/critical-fix`: Critical fixes for production

### Workflow

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Write clean, well-documented code
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**:
   ```bash
   pytest tests/
   black tmm/
   flake8 tmm/
   mypy tmm/
   ```

4. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add feature: descriptive commit message"
   ```

5. **Push and create a pull request**:
   ```bash
   git push origin feature/your-feature-name
   ```

## Code Standards

### Python Style

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) with some modifications:

- **Line length**: 88 characters (Black default)
- **Imports**: Use absolute imports, group by standard/third-party/local
- **Naming**: 
  - Functions and variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_CASE`
  - Private members: `_leading_underscore`

### Type Hints

All public functions and methods should include type hints:

```python
def calculate_reflectivity(
    energy_eV: Union[float, np.ndarray],
    angle_deg: float = 0.0,
    polarization: str = "s"
) -> Dict[str, np.ndarray]:
    """Calculate reflectivity spectrum."""
    ...
```

### Docstring Format

Use NumPy-style docstrings:

```python
def example_function(param1: int, param2: str) -> bool:
    """
    Brief description of the function.
    
    Longer description with more details about the function's
    purpose and behavior.
    
    Parameters
    ----------
    param1 : int
        Description of the first parameter
    param2 : str
        Description of the second parameter
        
    Returns
    -------
    bool
        Description of the return value
        
    Raises
    ------
    ValueError
        When invalid parameters are provided
        
    Examples
    --------
    >>> result = example_function(42, "test")
    >>> print(result)
    True
    """
    ...
```

### Error Handling

- Use specific exception types
- Provide informative error messages
- Include context about what went wrong and how to fix it

```python
if energy_eV <= 0:
    raise ValueError(
        f"Energy must be positive, got {energy_eV} eV. "
        "Please provide energy values greater than 0."
    )
```

## Testing

### Test Organization

Tests are organized in the `tests/` directory:

```
tests/
├── test_core.py              # Core TMM functionality
├── test_materials.py         # Material system tests
├── test_calculations.py      # Calculator tests
├── test_utils.py            # Utility function tests
└── test_integration.py      # Integration tests
```

### Writing Tests

- **Unit tests**: Test individual functions and methods
- **Integration tests**: Test component interactions
- **Property tests**: Test mathematical properties and invariants
- **Performance tests**: Benchmark critical functions

Example test:

```python
def test_material_creation():
    """Test basic material creation and properties."""
    energy = np.linspace(1.0, 3.0, 100)
    eps = 2.25 + 0.1j * np.ones_like(energy)
    
    material = Material("Test", energy, eps)
    
    assert material.name == "Test"
    assert len(material.energy_eV) == 100
    assert not material.is_anisotropic
    
    # Test interpolation
    eps_interp = material.get_dielectric_tensor(2.0)
    assert np.isclose(eps_interp[0], 2.25 + 0.1j)
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_materials.py

# Run with coverage
pytest --cov=tmm --cov-report=html

# Run performance tests
pytest tests/test_performance.py -v
```

## Documentation

### Code Documentation

- All public functions, classes, and methods must have docstrings
- Include type hints for all parameters and return values
- Provide examples in docstrings when helpful

### User Documentation

- **README.md**: Project overview and quick start
- **API documentation**: Generated from docstrings
- **Tutorials**: Step-by-step guides in Jupyter notebooks
- **Examples**: Practical use cases and demonstrations

### Building Documentation

```bash
# Install documentation dependencies
pip install -e ".[dev]"

# Build HTML documentation
cd docs/
make html

# View documentation
open _build/html/index.html
```

## Submitting Changes

### Pull Request Process

1. **Update documentation**: Ensure all changes are documented
2. **Add tests**: New features must include appropriate tests
3. **Update changelog**: Add entry to `CHANGELOG.md`
4. **Ensure CI passes**: All automated checks must pass
5. **Request review**: Assign appropriate reviewers

### Pull Request Template

```markdown
## Description
Brief description of changes made.

## Type of Change
- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes
- [ ] I have tested this change manually

## Checklist
- [ ] My code follows the style guidelines of this project
- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
```

### Commit Message Guidelines

Use conventional commit format:

```
type(scope): brief description

Longer description if needed

Fixes #123
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

## Community Guidelines

### Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please read our [Code of Conduct](CODE_OF_CONDUCT.md).

### Communication

- **GitHub Issues**: Bug reports, feature requests, and discussions
- **GitHub Discussions**: General questions and community support
- **Email**: For private or sensitive matters

### Getting Help

- Check the [documentation](https://tmm-optics.readthedocs.io)
- Search existing [issues](https://github.com/dingxu/tmm-optics/issues)
- Ask questions in [discussions](https://github.com/dingxu/tmm-optics/discussions)
- Join our community chat (link in README)

## Recognition

Contributors are recognized in:
- `AUTHORS.md` file
- Release notes
- Documentation acknowledgments

Thank you for contributing to the TMM package and helping advance optical simulation tools for the research community! 