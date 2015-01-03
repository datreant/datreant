==========================
Frequently Asked Questions
==========================

1. Why PyTables?

   `PyTables <https://github.com/PyTables/PyTables>`__ is a (fantastic)
   interface to the `hdf5 <http://www.hdfgroup.org/HDF5/>`__ data format.
   Although not itself a relational database, MDSynthesis uses PyTables for
   building and managing the persistent state files on disk for **Sim**,
   **Group**, and **Coordinator** objects. This was chosen over a traditional
   RDBS because we wanted MDSynthesis to be serverless, and SQLite was not
   ideal because its file locking mechanisms are known to be unreliable on a
   network file system (NFS).


