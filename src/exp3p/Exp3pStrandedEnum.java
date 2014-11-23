/*
 * Copyright 2014 Tomasz Konopka.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package exp3p;

/**
 * Strandedness of sequencing protocol
 * 
 * SAME will mean look for reads that are aligned onto the same strand as the gene def.
 * REVERSE will mean look for reads that are aligned on the opposite strand sas the gene def.
 * NONE will mean ignore the read strand, i.e. count all reads (default)
 * 
 * @author tkonopka
 */
public enum Exp3pStrandedEnum {
    SAME, REVERSE, NONE
}
