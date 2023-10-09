#pragma once

class OsiSolverInterface;

// @brief Compute the condition number in norm 1
double compute_condition_number_norm1(const OsiSolverInterface* const solver);

// @brief Compute the condition number in norm 2
double compute_condition_number_norm2(const OsiSolverInterface* const solver);

// @brief Compute the condition number in the given norm (1 or 2)
double compute_condition_number(const OsiSolverInterface* const solver, const int norm);
