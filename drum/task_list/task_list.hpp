#ifndef TASK_LIST_TASK_LIST_HPP_
#define TASK_LIST_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t
#include <string>       // std::string

// Athena++ headers
#include "../athena.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;
class FFTGravitySolverTaskList;
class TaskID;

// TODO(felker): these 4x declarations can be nested in TaskList if MGTaskList is derived

// constants = return codes for functions working on individual Tasks and TaskList
enum class TaskStatus {fail, success, next};
enum class TaskListStatus {running, stuck, complete, nothing_to_do};

//----------------------------------------------------------------------------------------
//! \class TaskID
//  \brief generalization of bit fields for Task IDs, status, and dependencies.

class TaskID {  // POD but not aggregate (there is a user-provided ctor)
 public:
  TaskID() = default;
  explicit TaskID(unsigned int id);
  void Clear();
  bool IsUnfinished(const TaskID& id) const;
  bool CheckDependencies(const TaskID& dep) const;
  void SetFinished(const TaskID& id);

  bool operator== (const TaskID& rhs) const;
  TaskID operator| (const TaskID& rhs) const;

 private:
  constexpr static int kNField_ = 2;
  std::uint64_t bitfld_[kNField_];

  friend class TaskList;
};


//----------------------------------------------------------------------------------------
//! \struct Task
//  \brief data and function pointer for an individual Task

struct Task { // aggregate and POD
  TaskID task_id;    // encodes task with bit positions in HydroIntegratorTaskNames
  TaskID dependency; // encodes dependencies to other tasks using " " " "
  TaskStatus (TaskList::*TaskFunc)(MeshBlock*, int);  // ptr to member function
  bool lb_time; // flag for automatic load balancing based on timing
};

//---------------------------------------------------------------------------------------
//! \struct TaskStates
//  \brief container for task states on a single MeshBlock

struct TaskStates { // aggregate and POD
  TaskID finished_tasks;
  int indx_first_task, num_tasks_left;
  void Reset(int ntasks) {
    indx_first_task = 0;
    num_tasks_left = ntasks;
    finished_tasks.Clear();
  }
};

//----------------------------------------------------------------------------------------
//! \class TaskList
//  \brief data and function definitions for task list base class

class TaskList {
 public:
  TaskList() : ntasks(0), nstages(0), task_list_{} {} // 2x direct + zero initialization
  // rule of five:
  virtual ~TaskList() = default;

  // data
  int ntasks;     // number of tasks in this list
  int nstages;    // number of times the tasklist is repeated per each full timestep

  // functions
  TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, int stage, TaskStates &ts);
  void DoTaskListOneStage(Mesh *pmesh, int stage);

 protected:
  // TODO(felker): rename to avoid confusion with class name
  Task task_list_[64*TaskID::kNField_];

 private:
  virtual void AddTask(const TaskID& id, const TaskID& dep) = 0;
  virtual void StartupTaskList(MeshBlock *pmb, int stage) = 0;
};

//----------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskList
//  \brief data and function definitions for TimeIntegratorTaskList derived class

class TimeIntegratorTaskList : public TaskList {
 public:
  TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm);

  //--------------------------------------------------------------------------------------
  //! \struct IntegratorWeight
  //  \brief weights used in time integrator tasks

  struct IntegratorWeight {
    // 2S or 3S* low-storage RK coefficients, Ketchenson (2010)
    Real delta; // low-storage coefficients to avoid double F() evaluation per substage
    Real gamma_1, gamma_2, gamma_3; // low-storage coeff for weighted ave of registers
    Real beta; // coeff. from bidiagonal Shu-Osher form Beta matrix, -1 diagonal terms
  };

  // data
  std::string integrator;
  Real cfl_limit; // dt stability limit for the particular time integrator + spatial order

  // functions
  TaskStatus ClearAllBoundary(MeshBlock *pmb, int stage);

  TaskStatus CalculateHydroFlux(MeshBlock *pmb, int stage);
  TaskStatus CalculateEMF(MeshBlock *pmb, int stage);
  TaskStatus CalculateRadiationFlux(MeshBlock *pmb, int stage);

  TaskStatus SendHydroFlux(MeshBlock *pmb, int stage);
  TaskStatus SendEMF(MeshBlock *pmb, int stage);

  TaskStatus ReceiveAndCorrectHydroFlux(MeshBlock *pmb, int stage);
  TaskStatus ReceiveAndCorrectEMF(MeshBlock *pmb, int stage);

  TaskStatus IntegrateHydro(MeshBlock *pmb, int stage);
  TaskStatus IntegrateField(MeshBlock *pmb, int stage);
  TaskStatus IntegrateChemistry(MeshBlock *pmb, int stage);

  TaskStatus AddSourceTermsHydro(MeshBlock *pmb, int stage);

  TaskStatus DiffuseHydro(MeshBlock *pmb, int stage);
  TaskStatus DiffuseField(MeshBlock *pmb, int stage);
  TaskStatus DiffuseScalars(MeshBlock *pmb, int stage);

  TaskStatus SendHydro(MeshBlock *pmb, int stage);
  TaskStatus SendField(MeshBlock *pmb, int stage);

  TaskStatus ReceiveHydro(MeshBlock *pmb, int stage);
  TaskStatus ReceiveField(MeshBlock *pmb, int stage);

  TaskStatus SetBoundariesHydro(MeshBlock *pmb, int stage);
  TaskStatus SetBoundariesField(MeshBlock *pmb, int stage);

  TaskStatus SendHydroShear(MeshBlock *pmb, int stage);
  TaskStatus ReceiveHydroShear(MeshBlock *pmb, int stage);
  TaskStatus SendFieldShear(MeshBlock *pmb, int stage);
  TaskStatus ReceiveFieldShear(MeshBlock *pmb, int stage);
  TaskStatus SendEMFShear(MeshBlock *pmb, int stage);
  TaskStatus ReceiveEMFShear(MeshBlock *pmb, int stage);
  TaskStatus RemapEMFShear(MeshBlock *pmb, int stage);

  TaskStatus Prolongation(MeshBlock *pmb, int stage);
  TaskStatus Primitives(MeshBlock *pmb, int stage);
  TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);
  TaskStatus UserWork(MeshBlock *pmb, int stage);
  TaskStatus NewBlockTimeStep(MeshBlock *pmb, int stage);
  TaskStatus CheckRefinement(MeshBlock *pmb, int stage);

  TaskStatus CalculateScalarFlux(MeshBlock *pmb, int stage);
  TaskStatus SendScalarFlux(MeshBlock *pmb, int stage);
  TaskStatus ReceiveScalarFlux(MeshBlock *pmb, int stage);
  TaskStatus IntegrateScalars(MeshBlock *pmb, int stage);
  TaskStatus SendScalars(MeshBlock *pmb, int stage);
  TaskStatus ReceiveScalars(MeshBlock *pmb, int stage);
  TaskStatus SetBoundariesScalars(MeshBlock *pmb, int stage);

  TaskStatus CalculateTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus SendTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus ReceiveTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus IntegrateTurbulence(MeshBlock *pmb, int stage);
  TaskStatus SendTurbulence(MeshBlock *pmb, int stage);
  TaskStatus ReceiveTurbulence(MeshBlock *pmb, int stage);
  TaskStatus SetBoundariesTurbulence(MeshBlock *pmb, int stage);

  // implicit tasks
  TaskStatus UpdateHydro(MeshBlock *pmb, int stage);

  // particle tasks
  TaskStatus IntegrateParticles(MeshBlock *pmb, int stage);
  TaskStatus MeshToParticles(MeshBlock *pmb, int stage);
  TaskStatus SendParticles(MeshBlock *pmb, int stage);
  TaskStatus ReceiveParticles(MeshBlock *pmb, int stage);
  TaskStatus AttachParticles(MeshBlock *pmb, int stage);
  TaskStatus ParticlesToMesh(MeshBlock *pmb, int stage);

 private:
  IntegratorWeight stage_wghts[MAX_NSTAGE];

  void AddTask(const TaskID& id, const TaskID& dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

//----------------------------------------------------------------------------------------
//! \class SuperTimeStepTaskList
//  \brief data and function definitions for SuperTimeStepTaskList derived class

class SuperTimeStepTaskList : public TaskList {
 public:
  SuperTimeStepTaskList(ParameterInput *pin, Mesh *pm, TimeIntegratorTaskList *ptlist);
  const Real sts_max_dt_ratio;

  // functions
  TaskStatus CalculateHydroFlux_STS(MeshBlock *pmb, int stage);
  TaskStatus CalculateEMF_STS(MeshBlock *pmb, int stage);
  TaskStatus CalculateScalarFlux_STS(MeshBlock *pmb, int stage);

  TaskStatus IntegrateHydro_STS(MeshBlock *pmb, int stage);
  TaskStatus IntegrateField_STS(MeshBlock *pmb, int stage);
  TaskStatus IntegrateScalars_STS(MeshBlock *pmb, int stage);

  TaskStatus PhysicalBoundary_STS(MeshBlock *pmb, int stage);

 private:
  // currently intiialized but unused. May use it for direct calls to TimeIntegrator fns:
  TimeIntegratorTaskList *ptlist_;
  void AddTask(const TaskID&, const TaskID& dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each hydro task.

// TODO(felker): uncomment the reserved TASK_NAMES once the features are merged to master
namespace HydroIntegratorTaskNames {

const TaskID NONE(0);
const TaskID CLEAR_ALLBND(1);

const TaskID CALC_HYDFLX(2);
const TaskID CALC_FLDFLX(3);
const TaskID CALC_RADFLX(4);
const TaskID CALC_CHMFLX(5);

const TaskID SEND_HYDFLX(6);
const TaskID SEND_FLDFLX(7);
// const TaskID SEND_RADFLX(8);
// const TaskID SEND_CHMFLX(9);

const TaskID RECV_HYDFLX(10);
const TaskID RECV_FLDFLX(11);
// const TaskID RECV_RADFLX(12);
// const TaskID RECV_CHMFLX(13);

const TaskID SRCTERM_HYD(14);
// const TaskID SRCTERM_FLD(15);
// const TaskID SRCTERM_RAD(16);
// const TaskID SRCTERM_CHM(17);

const TaskID INT_HYD(18);
const TaskID INT_FLD(19);
// const TaskID INT_RAD(20);
const TaskID INT_CHM(21);

const TaskID SEND_HYD(22);
const TaskID SEND_FLD(23);
// const TaskID SEND_RAD(24);
// const TaskID SEND_CHM(25);

const TaskID RECV_HYD(26);
const TaskID RECV_FLD(27);
// const TaskID RECV_RAD(28);
// const TaskID RECV_CHM(29);

const TaskID SETB_HYD(30);
const TaskID SETB_FLD(31);
// const TaskID SETB_RAD(32);
// const TaskID SETB_CHM(33);

const TaskID PROLONG(34);
const TaskID CONS2PRIM(35);
const TaskID PHY_BVAL(36);
const TaskID USERWORK(37);
const TaskID NEW_DT(38);
const TaskID FLAG_AMR(39);

const TaskID SEND_HYDSH(40);
const TaskID SEND_EMFSH(41);
const TaskID SEND_FLDSH(42);
const TaskID RECV_HYDSH(43);
const TaskID RECV_EMFSH(44);
const TaskID RECV_FLDSH(45);
const TaskID RMAP_EMFSH(46);

const TaskID DIFFUSE_HYD(47);
const TaskID DIFFUSE_FLD(48);

const TaskID CALC_SCLRFLX(49);
const TaskID SEND_SCLRFLX(50);
const TaskID RECV_SCLRFLX(51);
const TaskID INT_SCLR(52);
const TaskID SEND_SCLR(53);
const TaskID RECV_SCLR(54);
const TaskID SETB_SCLR(55);
const TaskID DIFFUSE_SCLR(56);

// const TaskID RECV_SCLRSH(57);
// const TaskID SEND_SCLRSH(58);

const TaskID CALC_TURBFLX(59);
const TaskID SEND_TURBFLX(60);
const TaskID RECV_TURBFLX(61);
const TaskID INT_TURB(62);
const TaskID SEND_TURB(63);
const TaskID RECV_TURB(64);
const TaskID SETB_TURB(65);

const TaskID UPDATE_HYD(69);
const TaskID INT_PART(70);
const TaskID MESH2PART(71);
const TaskID SEND_PART(72);
const TaskID RECV_PART(73);
const TaskID ATTACH_PART(74);
const TaskID PART2MESH(75);

}  // namespace HydroIntegratorTaskNames
#endif  // TASK_LIST_TASK_LIST_HPP_
