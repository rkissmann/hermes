#ifdef HERMES_HAVE_HDF5
#include "hermes/cosmicrays/Hdf5Reader.h"
#include <iostream>

Hdf5Reader::Hdf5Reader(std::string filename, int _suppress_verbosity)
{
	open_file(filename, _suppress_verbosity);
}

hid_t Hdf5Reader::open_file(std::string filename, int _suppress_verbosity)
{
	this->suppress_verbosity = _suppress_verbosity;

	hdf5file =  H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	h5group = H5Gopen2(hdf5file, "/Data", H5P_DEFAULT);

	return 0;

}

hid_t Hdf5Reader::close_file() {
	return H5Fclose(hdf5file);
}


template <typename T, typename>
hid_t Hdf5Reader::ReadGlobalAttribute(std::string AttrName, T &AttrData) {
  // Open Attribute
  hid_t info = H5Aopen(h5group, AttrName.c_str(), H5P_DEFAULT);
  // Read Attribute Data
  hid_t return_val = H5Aread(info, get_hdf5_data_type<T>(), &AttrData);
  H5Aclose(info);
  return return_val;
}


template <typename T, typename>
hid_t Hdf5Reader::ReadGlobalAttribute(std::string AttrName, std::vector<T> &AttrData) {
  // Open Attribute
  hid_t attr = H5Aopen(h5group, AttrName.c_str(), H5P_DEFAULT);

  // Get type of attribute
  hid_t atype  = H5Aget_type(attr);
  // Get attribute space
  hid_t aspace = H5Aget_space(attr);
  // Get rank of attribute
  int nEntries = H5Sget_simple_extent_npoints(aspace);

  AttrData.resize(nEntries);

  // Read Attribute Data
  hid_t return_val = H5Aread(attr, get_hdf5_data_type<T>(), &AttrData);
  H5Aclose(attr);
  return return_val;
}


template <typename T, typename>
hid_t Hdf5Reader::ReadAttributeFromDataset(std::string DataSetName,
		const std::string &AttrName, T &AttrData) {
  // Open the dataset
  hid_t dataset = H5Dopen2(h5group, DataSetName.c_str(), H5P_DEFAULT);

  hid_t attr = H5Aopen(dataset, AttrName.c_str(), H5P_DEFAULT);



  hid_t err = H5Aread(attr, get_hdf5_data_type<T>(), &AttrData);

  H5Aclose(attr);
  H5Dclose(dataset);

  return err;
}

#endif
