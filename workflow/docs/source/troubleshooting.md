# Troubleshooting


## Errors messages and solutions
### Error validating config file.

A problem in your `config.yaml` as been found. 

#### Forgotten a mandatory element 

Example :

```
ValidationError: 'allow_auto_import' is a required property
```
Solution : add `allow_auto_import: True` to `config.yaml`

#### Mispelled or Invalid element

Example :
`config.yaml`:
```yaml
...
alow_auto_import: True
...
```
```
ValidationError: Additional properties are not allowed ('alow_auto_import' was unexpected)
```
Solution : replace `alow_auto_import` with `allow_auto_import`
to `config.yaml`

#### Misstyped element

Example :

`config.yaml`:
```yaml
...
allow_auto_import : 45
...
```

```
ValidationError: 45 is not of type 'boolean'
Failed validating 'type' in schema['properties']['allow_auto_import']
```

Solution: Use True or False value
