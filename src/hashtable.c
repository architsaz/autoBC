#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hashtable.h"

// Hash function to convert a string key to an index
unsigned int hash(const char *key)
{
	unsigned int hash = 0;
	while (*key)
	{
		hash = (hash << 5) + (unsigned int)(*key++);
	}
	return hash % TABLE_SIZE;
}
// Custom strdup function if strdup is not available
char *my_strdup(const char *s)
{
	size_t len = strlen(s) + 1;
	char *dup = malloc(len);
	if (dup != NULL)
	{
		memcpy(dup, s, len);
	}
	return dup;
}
// Function to create a new entry
Entry *createEntry(const char *key, const char *value)
{
	Entry *entry = (Entry *)malloc(sizeof(Entry));
	entry->key = my_strdup(key);	 // Duplicate the string
	entry->value = my_strdup(value); // Duplicate the string
	entry->next = NULL;
	return entry;
}
// Function to insert a key-value pair into the hash table
void inserthash(HashTable *table, const char *key, const char *value)
{
	unsigned int index = hash(key);
	Entry *entry = table->entries[index];
	if (entry == NULL)
	{
		// No collision, simply insert the new entry
		table->entries[index] = createEntry(key, value);
	}
	else
	{
		// Handle collision by chaining
		while (entry->next != NULL && strcmp(entry->key, key) != 0)
		{
			entry = entry->next;
		}
		if (strcmp(entry->key, key) == 0)
		{
			// Update the value if the key already exists
			free(entry->value);
			entry->value = my_strdup(value);
		}
		else
		{
			// Add a new entry at the end of the chain
			entry->next = createEntry(key, value);
		}
	}
}
// Function to retrieve a value by key
char *gethash(HashTable *table, const char *key)
{
	unsigned int index = hash(key);
	Entry *entry = table->entries[index];
	while (entry != NULL)
	{
		if (strcmp(entry->key, key) == 0)
		{
			return entry->value;
		}
		entry = entry->next;
	}
	return NULL; // Key not found
}
// Function to free the memory allocated for the hash table
void freeTable(HashTable *table)
{
	for (int i = 0; i < TABLE_SIZE; i++)
	{
		Entry *entry = table->entries[i];
		while (entry != NULL)
		{
			Entry *temp = entry;
			entry = entry->next;
			free(temp->key);
			free(temp->value);
			free(temp);
		}
	}
}
