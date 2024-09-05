include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | verify_data_structure.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset": checks["exit_code"] == 0 ? state.input : null,
        ]
      }
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset != null
    }

    | select_reference.run(
      fromState: [
        input: "dataset",
        num_features: "num_reference_genes",
        coord_type_proc: "coord_type_proc"
      ],
      toState: [dataset: "output"]
    )

    | simulate_svg.run(
      fromState: [
        input: "dataset",
        gp_k: "gp_k_sim",
        select_top_variable_genes: "select_top_variable_genes_sim"
      ],
      toState: [
        dataset_simulated: "output"
      ]
    )

    | log_cp.run(
      fromState: [
        input: "dataset_simulated",
      ],
      toState: [
        dataset_simulated_normalized: "output"
      ],
      args: [n_cp: -1]
    )

    | split_dataset.run(
      fromState: [
        input: "dataset_simulated_normalized"
      ],
      toState: [
        output_dataset: "output_dataset",
        output_solution: "output_solution" 
      ]
    )

    // only output the files for which an output file was specified
    | setState(["output_dataset", "output_solution", "dataset_simulated_normalized"])

  emit:
  output_ch
}


// temp fix for rename_keys typo

def findStatesTemp(Map params, Map config) {
  def auto_config = deepClone(config)
  def auto_params = deepClone(params)

  auto_config = auto_config.clone()
  // override arguments
  auto_config.argument_groups = []
  auto_config.arguments = [
    [
      type: "string",
      name: "--id",
      description: "A dummy identifier",
      required: false
    ],
    [
      type: "file",
      name: "--input_states",
      example: "/path/to/input/directory/**/state.yaml",
      description: "Path to input directory containing the datasets to be integrated.",
      required: true,
      multiple: true,
      multiple_sep: ";"
    ],
    [
      type: "string",
      name: "--filter",
      example: "foo/.*/state.yaml",
      description: "Regex to filter state files by path.",
      required: false
    ],
    // to do: make this a yaml blob?
    [
      type: "string",
      name: "--rename_keys",
      example: ["newKey1:oldKey1", "newKey2:oldKey2"],
      description: "Rename keys in the detected input files. This is useful if the input files do not match the set of input arguments of the workflow.",
      required: false,
      multiple: true,
      multiple_sep: ";"
    ],
    [
      type: "string",
      name: "--settings",
      example: '{"output_dataset": "dataset.h5ad", "k": 10}',
      description: "Global arguments as a JSON glob to be passed to all components.",
      required: false
    ]
  ]
  if (!(auto_params.containsKey("id"))) {
    auto_params["id"] = "auto"
  }

  // run auto config through processConfig once more
  auto_config = processConfig(auto_config)

  workflow findStatesTempWf {
    helpMessage(auto_config)

    output_ch = 
      channelFromParams(auto_params, auto_config)
        | flatMap { autoId, args ->

          def globalSettings = args.settings ? readYamlBlob(args.settings) : [:]

          // look for state files in input dir
          def stateFiles = args.input_states

          // filter state files by regex
          if (args.filter) {
            stateFiles = stateFiles.findAll{ stateFile ->
              def stateFileStr = stateFile.toString()
              def matcher = stateFileStr =~ args.filter
              matcher.matches()}
          }

          // read in states
          def states = stateFiles.collect { stateFile ->
            def state_ = readTaggedYaml(stateFile)
            [state_.id, state_]
          }

          // construct renameMap
          if (args.rename_keys) {
            def renameMap = args.rename_keys.collectEntries{renameString ->
              def split = renameString.split(":")
              assert split.size() == 2: "Argument 'rename_keys' should be of the form 'newKey:oldKey;newKey:oldKey'"
              split
            }

            // rename keys in state, only let states through which have all keys
            // also add global settings
            states = states.collectMany{id, state ->
              def newState = [:]

              for (key in renameMap.keySet()) {
                def origKey = renameMap[key]
                if (!(state.containsKey(origKey))) {
                  return []
                }
                newState[key] = state[origKey]
              }

              [[id, globalSettings + newState]]
            }
          }

          states
        }
    emit:
    output_ch
  }

  return findStatesTempWf
}