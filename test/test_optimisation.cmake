# Integration test: 2pF Fermi parameter optimisation for Pb-206
# Expected result: rms_radius = 5.4853 fm (Anderson et al., Phys. Rev. 187, 1565 (1969))
# Tolerance: ±0.05 fm

execute_process(
  COMMAND "${MUDIRAC_BIN}" pb206_lm.in pb206_exp.in
  WORKING_DIRECTORY "${DATA_DIR}"
  RESULT_VARIABLE exit_code
)
if(NOT exit_code EQUAL 0)
  message(FATAL_ERROR "mudirac exited with code ${exit_code}")
endif()

# Output file has three lines: comment, header, data
# Data columns: fermi_c  fermi_t  rms_radius  theta  mean_chi_sq  time
file(STRINGS "${DATA_DIR}/pb206_lm.fermi_parameters.out" file_lines)
list(GET file_lines 2 data_line)
string(REGEX MATCHALL "[0-9]+\\.[0-9]+" values "${data_line}")
list(GET values 2 rms_radius)

message(STATUS "Pb-206 optimisation: rms_radius = ${rms_radius} fm (expected ~5.4853 fm)")

if(rms_radius LESS 5.45 OR rms_radius GREATER 5.52)
  message(FATAL_ERROR "rms_radius ${rms_radius} fm outside tolerance [5.45, 5.52] fm")
endif()
