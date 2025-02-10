#ifndef HASHTABLE_H
#define HASHTABLE_H

#define TABLE_SIZE 500
// Define the structure for a hash table entry
typedef struct Entry
{
    char *key;
    char *value;
    struct Entry *next; // For handling collisions using chaining
} Entry;
// Define the structure for the hash table
typedef struct
{
    Entry *entries[TABLE_SIZE];
} HashTable;
void inserthash(HashTable *table, const char *key, const char *value);
char *gethash(HashTable *table, const char *key);
void freeTable(HashTable *table);
#endif