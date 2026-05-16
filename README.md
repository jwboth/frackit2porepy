# frackit2porepy
Script to generate fracture networks using Frackit, and exported to PorePy compatible format.

To run the fracture generator:
```bash
docker build -t frackit2porepy:latest .
docker run --rm -v ".\shared:/frackit/shared" frackit2porepy:latest
```

The generator can be configured in the shared config file `shared/config.toml`.
Follow the instructions in the comments of the configuration file.