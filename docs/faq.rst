==========================
Frequently Asked Questions
==========================

1. Why PyTables?

   `PyTables <https://github.com/PyTables/PyTables>`__ is a (fantastic)
   interface to the `hdf5 <http://www.hdfgroup.org/HDF5/>`__ data format.
   Although not itself a relational database, datreant uses PyTables for
   building and managing the persistent state files on disk for Treants. This
   was chosen over a traditional RDBS because we wanted datreant to be
   serverless. An alternative format for state files as SQLite databases
   may be developed in the future.


