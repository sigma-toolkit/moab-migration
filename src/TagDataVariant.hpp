#ifndef MOAB_TAG_DATA_VARIANT_HPP
#define MOAB_TAG_DATA_VARIANT_HPP

#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include <vector>
#include <memory>
#include <type_traits>
#include <functional>

namespace moab {

/**
 * @brief Type-safe variant for MOAB tag data with automatic type conversion
 *
 * This class implements a visitor pattern to handle different numeric types
 * commonly used in MOAB tags, providing automatic conversion between types
 * during get/set operations.
 */
class TagDataVariant {
public:
    enum DataTypeEnum {
        TYPE_INT,
        TYPE_UNSIGNED_INT,
        TYPE_LONG,
        TYPE_UNSIGNED_LONG,
        TYPE_UNSIGNED_LONG_LONG,
        TYPE_UNKNOWN
    };

    /**
     * @brief Visit pattern - applies visitor to the underlying typed vector
     */
    template<typename Visitor>
    auto visit(Visitor&& visitor) -> decltype(auto) {
        switch(type_) {
            case TYPE_INT:
                return visitor(*static_cast<std::vector<int>*>(data_.get()));
            case TYPE_UNSIGNED_INT:
                return visitor(*static_cast<std::vector<unsigned int>*>(data_.get()));
            case TYPE_LONG:
                return visitor(*static_cast<std::vector<long>*>(data_.get()));
            case TYPE_UNSIGNED_LONG:
                return visitor(*static_cast<std::vector<unsigned long>*>(data_.get()));
            case TYPE_UNSIGNED_LONG_LONG:
                return visitor(*static_cast<std::vector<unsigned long long>*>(data_.get()));
            default:
                throw std::runtime_error("Invalid data type for visit");
        }
    }

    /**
     * @brief Visit pattern - applies visitor to the underlying typed vector (const version)
     */
    template<typename Visitor>
    auto visit(Visitor&& visitor) const -> decltype(auto) {
        switch(type_) {
            case TYPE_INT:
                return visitor(*static_cast<const std::vector<int>*>(data_.get()));
            case TYPE_UNSIGNED_INT:
                return visitor(*static_cast<const std::vector<unsigned int>*>(data_.get()));
            case TYPE_LONG:
                return visitor(*static_cast<const std::vector<long>*>(data_.get()));
            case TYPE_UNSIGNED_LONG:
                return visitor(*static_cast<const std::vector<unsigned long>*>(data_.get()));
            case TYPE_UNSIGNED_LONG_LONG:
                return visitor(*static_cast<const std::vector<unsigned long long>*>(data_.get()));
            default:
                throw std::runtime_error("Invalid data type for visit");
        }
    }

private:
    DataTypeEnum type_;
    std::unique_ptr<void, std::function<void(void*)>> data_;
    size_t size_;

    // Type trait mapping from MOAB DataType to our enum
    static DataTypeEnum moab_type_to_enum(DataType dt) {
        switch(dt) {
            case MB_TYPE_INTEGER: return TYPE_INT;
            case MB_TYPE_UNSIGNED_INTEGER: return TYPE_UNSIGNED_INT;
            case MB_TYPE_LONG: return TYPE_LONG;
            case MB_TYPE_UNSIGNED_LONG: return TYPE_UNSIGNED_LONG;
            case MB_TYPE_UNSIGNED_LONG_LONG: return TYPE_UNSIGNED_LONG_LONG;
            default: return TYPE_UNKNOWN;
        }
    }

    // Type trait mapping from C++ types to our enum (C++11 compatible)
    template<typename T>
    static DataTypeEnum cpp_type_to_enum() {
        if (std::is_same<T, int>::value) return TYPE_INT;
        else if (std::is_same<T, unsigned int>::value) return TYPE_UNSIGNED_INT;
        else if (std::is_same<T, long>::value) return TYPE_LONG;
        else if (std::is_same<T, unsigned long>::value) return TYPE_UNSIGNED_LONG;
        else if (std::is_same<T, unsigned long long>::value) return TYPE_UNSIGNED_LONG_LONG;
        else return TYPE_UNKNOWN;
    }

    // Create typed vector and return void pointer with custom deleter
    template<typename T>
    std::unique_ptr<void, std::function<void(void*)>> create_vector(size_t size) {
        auto* vec = new std::vector<T>(size);
        return std::unique_ptr<void, std::function<void(void*)>>(
            vec,
            [](void* ptr) { delete static_cast<std::vector<T>*>(ptr); }
        );
    }

public:
    /**
     * @brief Construct from MOAB tag information
     */
    TagDataVariant(Interface* mb, Tag tag, size_t num_entities) {
        DataType moab_dtype;
        int tag_length;
        mb->tag_get_data_type(tag, moab_dtype);
        mb->tag_get_length(tag, tag_length);

        type_ = moab_type_to_enum(moab_dtype);
        size_ = num_entities * tag_length;

        // Get the actual tag size in bytes
        int tag_bytes;
        mb->tag_get_bytes(tag, tag_bytes);

        // Allocate buffer with proper size based on the target type
        // We need to ensure the buffer is large enough to hold the data in the target type
        size_t elements_needed = size_;

        // If we're reading into a larger type, we need to adjust the number of elements
        size_t type_size = 0;
        switch(type_) {
            case TYPE_INT:              type_size = sizeof(int); break;
            case TYPE_UNSIGNED_INT:     type_size = sizeof(unsigned int); break;
            case TYPE_LONG:             type_size = sizeof(long); break;
            case TYPE_UNSIGNED_LONG:    type_size = sizeof(unsigned long); break;
            case TYPE_UNSIGNED_LONG_LONG: type_size = sizeof(unsigned long long); break;
            default:
                throw std::runtime_error("Unsupported tag data type");
        }

        // Ensure we have enough space for the target type
        elements_needed = (elements_needed * tag_bytes + type_size - 1) / type_size;

        // Allocate buffer with proper size
        switch(type_) {
            case TYPE_INT:              data_ = create_vector<int>(elements_needed); break;
            case TYPE_UNSIGNED_INT:     data_ = create_vector<unsigned int>(elements_needed); break;
            case TYPE_LONG:             data_ = create_vector<long>(elements_needed); break;
            case TYPE_UNSIGNED_LONG:    data_ = create_vector<unsigned long>(elements_needed); break;
            case TYPE_UNSIGNED_LONG_LONG: data_ = create_vector<unsigned long long>(elements_needed); break;
            default:
                throw std::runtime_error("Unsupported tag data type");
        }

        // Initialize the buffer to zero to prevent any uninitialized reads
        visit([](auto& vec) { std::fill(vec.begin(), vec.end(), typename std::remove_reference<decltype(vec)>::type::value_type(0)); });
    }

    /**
     * @brief Construct with explicit type
     */
    template<typename T>
    TagDataVariant(size_t size) : type_(cpp_type_to_enum<T>()), size_(size) {
        data_ = create_vector<T>(size);
    }

    /**
     * @brief Get raw data pointer for MOAB interface
     * @note The returned pointer is valid as long as this object exists
     * and no non-const methods are called on it.
     */
    void* data() {
        switch(type_) {
            case TYPE_INT:              return static_cast<std::vector<int>*>(data_.get())->data();
            case TYPE_UNSIGNED_INT:     return static_cast<std::vector<unsigned int>*>(data_.get())->data();
            case TYPE_LONG:             return static_cast<std::vector<long>*>(data_.get())->data();
            case TYPE_UNSIGNED_LONG:    return static_cast<std::vector<unsigned long>*>(data_.get())->data();
            case TYPE_UNSIGNED_LONG_LONG: return static_cast<std::vector<unsigned long long>*>(data_.get())->data();
            default: return nullptr;
        }
    }

    /**
     * @brief Get raw data pointer for MOAB interface (const version)
     * @note The returned pointer is valid as long as this object exists
     * and no non-const methods are called on it.
     */
    const void* data() const {
        switch(type_) {
            case TYPE_INT:              return static_cast<const std::vector<int>*>(data_.get())->data();
            case TYPE_UNSIGNED_INT:     return static_cast<const std::vector<unsigned int>*>(data_.get())->data();
            case TYPE_LONG:             return static_cast<const std::vector<long>*>(data_.get())->data();
            case TYPE_UNSIGNED_LONG:    return static_cast<const std::vector<unsigned long>*>(data_.get())->data();
            case TYPE_UNSIGNED_LONG_LONG: return static_cast<const std::vector<unsigned long long>*>(data_.get())->data();
            default: return nullptr;
        }
    }

    /**
     * @brief Get the number of bytes in the data buffer
     */
    size_t data_size() const {
        switch(type_) {
            case TYPE_INT:              return static_cast<const std::vector<int>*>(data_.get())->size() * sizeof(int);
            case TYPE_UNSIGNED_INT:     return static_cast<const std::vector<unsigned int>*>(data_.get())->size() * sizeof(unsigned int);
            case TYPE_LONG:             return static_cast<const std::vector<long>*>(data_.get())->size() * sizeof(long);
            case TYPE_UNSIGNED_LONG:    return static_cast<const std::vector<unsigned long>*>(data_.get())->size() * sizeof(unsigned long);
            case TYPE_UNSIGNED_LONG_LONG: return static_cast<const std::vector<unsigned long long>*>(data_.get())->size() * sizeof(unsigned long long);
            default: return 0;
        }
    }

    size_t size() const { return size_; }
    DataTypeEnum type() const { return type_; }

    // Copy data to target vector with type conversion
    template<typename TargetType>
    void copy_to(std::vector<TargetType>& target) const {
        target.resize(size_);

        auto copy_impl = [&](const auto* vec) {
            using SourceType = typename std::remove_pointer_t<decltype(vec)>::value_type;
            size_t elements_to_copy = std::min({(size_t)size_, vec->size(), target.size()});
            if (elements_to_copy > 0) {
                std::transform(vec->begin(),
                             vec->begin() + elements_to_copy,
                             target.begin(),
                             [](const SourceType& val) {
                                 return static_cast<TargetType>(val);
                             });
            }
            if (elements_to_copy < target.size()) {
                std::fill(target.begin() + elements_to_copy, target.end(), TargetType(0));
            }
        };

        switch(type_) {
            case TYPE_INT:              copy_impl(static_cast<const std::vector<int>*>(data_.get())); break;
            case TYPE_UNSIGNED_INT:     copy_impl(static_cast<const std::vector<unsigned int>*>(data_.get())); break;
            case TYPE_LONG:             copy_impl(static_cast<const std::vector<long>*>(data_.get())); break;
            case TYPE_UNSIGNED_LONG:    copy_impl(static_cast<const std::vector<unsigned long>*>(data_.get())); break;
            case TYPE_UNSIGNED_LONG_LONG: copy_impl(static_cast<const std::vector<unsigned long long>*>(data_.get())); break;
            default: throw std::runtime_error("Unsupported data type in copy_to");
        }
    }

    // Copy data from source vector with type conversion
    template<typename SourceType>
    void copy_from(const std::vector<SourceType>& source) {
        if (source.size() != size_) {
            throw std::runtime_error("Size mismatch in copy_from: expected " +
                                  std::to_string(size_) + " elements, got " +
                                  std::to_string(source.size()));
        }

        auto copy_impl = [&](auto* vec) {
            using TargetType = typename std::remove_pointer_t<decltype(vec)>::value_type;
            size_t elements_to_copy = std::min({(size_t)size_, vec->size(), source.size()});

            if (elements_to_copy > 0) {
                std::transform(source.begin(),
                             source.begin() + elements_to_copy,
                             vec->begin(),
                             [](const SourceType& val) {
                                 return static_cast<TargetType>(val);
                             });
            }

            if (elements_to_copy < vec->size()) {
                std::fill(vec->begin() + elements_to_copy, vec->end(), TargetType(0));
            }
        };

        switch(type_) {
            case TYPE_INT:              copy_impl(static_cast<std::vector<int>*>(data_.get())); break;
            case TYPE_UNSIGNED_INT:     copy_impl(static_cast<std::vector<unsigned int>*>(data_.get())); break;
            case TYPE_LONG:             copy_impl(static_cast<std::vector<long>*>(data_.get())); break;
            case TYPE_UNSIGNED_LONG:    copy_impl(static_cast<std::vector<unsigned long>*>(data_.get())); break;
            case TYPE_UNSIGNED_LONG_LONG: copy_impl(static_cast<std::vector<unsigned long long>*>(data_.get())); break;
            default: throw std::runtime_error("Unsupported data type in copy_from");
        }
    }
};

/**
 * @brief High-level interface for type-safe tag operations
 */
class TypeSafeTagInterface {
private:
    Interface* mb_;
    Tag tag_;
    DataType dtype_;
    int length_;

public:
    TypeSafeTagInterface(Interface* mb, Tag tag) : mb_(mb), tag_(tag) {
        mb_->tag_get_data_type(tag_, dtype_);
        mb_->tag_get_length(tag_, length_);
    }

    TypeSafeTagInterface(Interface* mb, const std::string& tag_name) : mb_(mb) {
        ErrorCode rval = mb_->tag_get_handle(tag_name.c_str(), tag_);
        if (MB_SUCCESS != rval) {
            throw std::runtime_error("Tag '" + tag_name + "' not found");
        }
        mb_->tag_get_data_type(tag_, dtype_);
        mb_->tag_get_length(tag_, length_);
    }

    /**
     * @brief Get tag data with automatic type conversion
     */
    template<typename TargetType>
    ErrorCode get_data(const Range& entities, std::vector<TargetType>& values) {
        if (entities.empty()) {
            values.clear();
            return MB_SUCCESS;
        }

        // Create a variant with proper type and size
        TagDataVariant variant(mb_, tag_, entities.size());

        // Get data from MOAB into variant
        ErrorCode rval = mb_->tag_get_data(tag_, entities, variant.data());
        if (MB_SUCCESS != rval) {
            values.clear();
            return rval;
        }

        try {
            // Convert to target type
            variant.copy_to(values);
            return MB_SUCCESS;
        } catch (const std::exception& e) {
            values.clear();
            return MB_FAILURE;
        }
    }

    /**
     * @brief Get tag data with automatic type conversion
     */
    template<typename TargetType>
    ErrorCode get_data(const std::vector<EntityHandle>& entities, std::vector<TargetType>& values) {
        if (entities.empty()) {
            values.clear();
            return MB_SUCCESS;
        }

        // Create a variant with proper type and size
        TagDataVariant variant(mb_, tag_, entities.size());

        // Get data from MOAB into variant
        ErrorCode rval = mb_->tag_get_data(tag_, entities.data(), entities.size(), variant.data());
        if (MB_SUCCESS != rval) {
            values.clear();
            return rval;
        }

        try {
            // Convert to target type
            variant.copy_to(values);
            return MB_SUCCESS;
        } catch (const std::exception& e) {
            values.clear();
            return MB_FAILURE;
        }
    }

    /**
     * @brief Get tag data for single entity
     */
    template<typename TargetType>
    ErrorCode get_data(EntityHandle entity, std::vector<TargetType>& values) {
        Range single_entity;
        single_entity.insert(entity);
        return get_data(single_entity, values);
    }

    /**
     * @brief Get tag data for single entity
     */
    template<typename TargetType>
    ErrorCode get_data(EntityHandle entity, TargetType& value) {
        Range single_entity;
        single_entity.insert(entity);
        std::vector<TargetType> values(1);
        MB_CHK_SET_ERR(get_data(single_entity, values), "can't get tag data");
        value = values[0];
        return MB_SUCCESS;
    }

    /**
     * @brief Set tag data with automatic type conversion
     * @param entities The entities to set tag data for
     * @param values Input vector containing the values to set
     * @return ErrorCode indicating success or failure
     * @note The values vector should contain tag_length values for each entity,
     *       stored contiguously (entity1_val1, entity1_val2, ..., entity2_val1, ...).
     *       The size of values must be exactly entities.size() * tag_length.
     */
    template<typename SourceType>
    ErrorCode set_data(const Range& entities, const std::vector<SourceType>& values) {
        if (entities.empty()) {
            return MB_SUCCESS;  // Nothing to do for empty range
        }

        if (values.size() != entities.size() * length_) {
            return MB_FAILURE;  // Size mismatch
        }

        try {
            // Create variant and copy data
            TagDataVariant variant(mb_, tag_, entities.size());
            variant.copy_from(values);

            // Set data in MOAB
            return mb_->tag_set_data(tag_, entities, variant.data());
        } catch (const std::exception& e) {
            return MB_FAILURE;
        }
    }

    /**
     * @brief Set tag data with automatic type conversion
     * @param entities The entities to set tag data for
     * @param values Input vector containing the values to set
     * @return ErrorCode indicating success or failure
     * @note The values vector should contain tag_length values for each entity,
     *       stored contiguously (entity1_val1, entity1_val2, ..., entity2_val1, ...).
     *       The size of values must be exactly entities.size() * tag_length.
     */
    template<typename SourceType>
    ErrorCode set_data(const std::vector<EntityHandle>& entities, const std::vector<SourceType>& values) {
        if (entities.empty()) {
            return MB_SUCCESS;  // Nothing to do for empty range
        }

        if (values.size() != entities.size() * length_) {
            return MB_FAILURE;  // Size mismatch
        }

        try {
            // Create variant and copy data
            TagDataVariant variant(mb_, tag_, entities.size());
            variant.copy_from(values);

            // Set data in MOAB
            return mb_->tag_set_data(tag_, entities.data(), entities.size(), variant.data());
        } catch (const std::exception& e) {
            return MB_FAILURE;
        }
    }

    /**
     * @brief Set tag data for single entity
     */
    template<typename SourceType>
    ErrorCode set_data(EntityHandle entity, const std::vector<SourceType>& values) {
        Range single_entity;
        single_entity.insert(entity);
        return set_data(single_entity, values);
    }

    /**
     * @brief Set tag data for single entity
     */
    template<typename SourceType>
    ErrorCode set_data(EntityHandle entity, const SourceType& value) {
        Range single_entity;
        single_entity.insert(entity);
        std::vector<SourceType> values(1, value);
        MB_CHK_SET_ERR(set_data(single_entity, values), "can't set tag data");
        return MB_SUCCESS;
    }

    /**
     * @brief Get tag information
     */
    DataType get_data_type() const { return dtype_; }
    int get_length() const { return length_; }
    Tag get_tag() const { return tag_; }

    /**
     * @brief Create a variant buffer for this tag
     */
    TagDataVariant create_buffer(size_t num_entities) {
        return TagDataVariant(mb_, tag_, num_entities);
    }
};

/**
 * @brief Factory function for easy creation
 */
inline TypeSafeTagInterface make_tag_interface(Interface* mb, const std::string& tag_name) {
    return TypeSafeTagInterface(mb, tag_name);
}

inline TypeSafeTagInterface make_tag_interface(Interface* mb, Tag tag) {
    return TypeSafeTagInterface(mb, tag);
}

} // namespace moab

#endif // MOAB_TAG_DATA_VARIANT_HPP
