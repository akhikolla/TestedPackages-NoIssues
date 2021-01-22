
#ifndef MACROS_H
#define MACROS_H




/**
 * @author	Wilson Toussile
 * @brief Macros
 */

#define MyError(message) Rprintf("\n >>>> Error : %s in %s\n", message, __PRETTY_FUNCTION__)
#define MyDebug(message) Rprintf("\n >>> Debug %s in %s\n", message, __PRETTY_FUNCTION__)
#define MyWarning(message) Rprintf("\n >> Wraning : %s in %s\n", message, __PRETTY_FUNCTION__)

#define MyWarningVar(varName, message) Rprintf("\n >> Warning : Argument %s in %s : %s\n", varName, __PRETTY_FUNCTION__, message)

#define ErrorMemoryAlloc() Rprintf("\n >>> Memory allocation problem in %s\n", __PRETTY_FUNCTION__);
#define ErrorIndexOutOfRange() Rprintf("\n >>> Index out of range in %s\n", __PRETTY_FUNCTION__);
#define ErrorOpeningFile() Rprintf("\n >>> Unable to open file %s\n", __PRETTY_FUNCTION__);
#define ErrorReadingStream() Rprintf("\n >>> Error reading OIStream in %s", __PRETTY_FUNCTION__);
#define ErrorInvalideParameter() Rprintf("\n >>> Invalide parameter in %s", __PRETTY_FUNCTION__);
#define ErrorInvalideDimensions() MyError("\n >>> Invalide dimensions in ", __PRETTY_FUNCTION__);
#define ErrorObjectNotAllocated() MyError("\n >>> Object not allocated in ", __PRETTY_FUNCTION__);


#define MyDebugMessage() Rprintf("Debug in  %s", __PRETTY_FUNCTION__);

/** Variable name */
#define MyVarName(x) #x

/** -----------------------------------------------------------*/

#endif /* MACROS_H*/

